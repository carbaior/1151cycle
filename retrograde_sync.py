#!/usr/bin/env python3
"""
retrograde_sync.py
==================
Quantifies the synchronisation of planetary retrograde motions between
the reference series (starting at t0) and the candidate series displaced
by T* = 420,403 days.

For each outer planet (Mars, Jupiter, Saturn) and inner planet (Mercury,
Venus), identifies retrograde episodes in both series and measures the
temporal shift between corresponding episodes.

Retrograde definition (geocentric):
  A planet is in retrograde when its geocentric ecliptic longitude
  decreases from one day to the next.

Output:
  retrograde_sync.csv   — per-episode comparison table
  retrograde_stats.csv  — summary statistics per planet

Requires: helio6_de441_cache.pkl  (from helio1151.py)
          Skyfield + DE441 ephemeris for geocentric positions

Author: Carlos Baiget Orts — asinfreedom@gmail.com
"""

import pickle, csv, os
import numpy as np
from skyfield.api import load, GREGORIAN_START
from skyfield.framelib import ecliptic_J2000_frame

# ── Configuration ─────────────────────────────────────────────────────────────

EPHEMERIS   = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp'
GEO_CACHE   = 'geo_cache.pkl'          # geocentric cache (built here if absent)

T_STAR      = 420403                   # days
SERIES_DAYS = 36525                    # 100 Julian years

REF_YEAR, REF_MONTH, REF_DAY = 0, 6, 15

PLANET_NAMES   = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Neptune']
PLANET_TARGETS = [
    'mercury barycenter',
    'venus barycenter',
    'mars barycenter',
    'jupiter barycenter',
    'saturn barycenter',
    'neptune barycenter',
]

# Minimum retrograde duration to count as a genuine episode (days)
MIN_RETRO_DAYS = {
    'Mercury': 18,
    'Venus':   30,
    'Mars':    50,
    'Jupiter': 80,
    'Saturn':  90,
    'Neptune': 120,
}

# ── Build geocentric position cache ──────────────────────────────────────────

def build_geo_cache(ts, eph, jd_ref, series_days, t_star):
    """
    Compute daily geocentric ecliptic longitudes for the 5 planets
    over two series: reference [jd_ref, jd_ref+N] and
    candidate [jd_ref-T*, jd_ref-T*+N].
    Returns dict: planet -> (lon_ref array, lon_cand array), shape (N,)
    """
    earth   = eph['earth']
    targets = [eph[t] for t in PLANET_TARGETS]

    n = series_days
    result = {}

    for pi, (name, tgt) in enumerate(zip(PLANET_NAMES, targets)):
        print(f"  Computing geocentric positions: {name}...", end='\r')
        lon_ref  = np.zeros(n, dtype=np.float64)
        lon_cand = np.zeros(n, dtype=np.float64)

        batch = 500
        for start in range(0, n, batch):
            end = min(start + batch, n)

            # Reference series
            jds_r = np.arange(jd_ref + start, jd_ref + end, dtype=float)
            t_r   = ts.tt_jd(jds_r)
            ast_r = earth.at(t_r).observe(tgt).apparent()
            _, lon_r, _ = ast_r.frame_latlon(ecliptic_J2000_frame)
            lon_ref[start:end] = lon_r.degrees

            # Candidate series (T* earlier)
            jds_c = np.arange(jd_ref - t_star + start,
                              jd_ref - t_star + end, dtype=float)
            t_c   = ts.tt_jd(jds_c)
            ast_c = earth.at(t_c).observe(tgt).apparent()
            _, lon_c, _ = ast_c.frame_latlon(ecliptic_J2000_frame)
            lon_cand[start:end] = lon_c.degrees

        result[name] = (lon_ref, lon_cand)
        print(f"  {name}: done.             ")

    return result

def load_or_build_geo(ts, eph, jd_ref, series_days, t_star):
    if os.path.exists(GEO_CACHE):
        print(f"Loading geocentric cache from {GEO_CACHE}...")
        with open(GEO_CACHE, 'rb') as f:
            data = pickle.load(f)
        if (data.get('jd_ref') == jd_ref and
                data.get('series_days') == series_days and
                data.get('t_star') == t_star):
            print("Cache OK.")
            return data['positions']
        print("Cache parameters changed — rebuilding.")

    print("Building geocentric position cache...")
    positions = build_geo_cache(ts, eph, jd_ref, series_days, t_star)
    with open(GEO_CACHE, 'wb') as f:
        pickle.dump({
            'jd_ref': jd_ref, 'series_days': series_days,
            't_star': t_star, 'positions': positions
        }, f)
    print(f"Cache saved to {GEO_CACHE}")
    return positions

# ── Retrograde detection ──────────────────────────────────────────────────────

def circular_diff(a, b):
    """Signed difference a-b on circle, result in (-180, 180]."""
    d = (a - b) % 360.0
    return np.where(d > 180.0, d - 360.0, d)

def find_retrogrades(lon_series, min_days):
    """
    Find retrograde episodes in a longitude series.
    Returns list of (start_day, end_day, duration_days).
    """
    n = len(lon_series)
    # Daily motion: positive = prograde, negative = retrograde
    motion = circular_diff(lon_series[1:], lon_series[:-1])

    episodes = []
    in_retro  = False
    start_day = 0

    for i in range(len(motion)):
        if not in_retro and motion[i] < 0:
            in_retro  = True
            start_day = i
        elif in_retro and motion[i] >= 0:
            dur = i - start_day
            if dur >= min_days:
                # peak = day of maximum retrograde speed
                peak = start_day + np.argmin(motion[start_day:i])
                episodes.append({
                    'start': start_day,
                    'end':   i,
                    'dur':   dur,
                    'peak':  int(peak),
                })
            in_retro = False

    # Handle episode that reaches end of series
    if in_retro:
        dur = len(motion) - start_day
        if dur >= min_days:
            peak = start_day + np.argmin(motion[start_day:])
            episodes.append({
                'start': start_day,
                'end':   len(motion),
                'dur':   dur,
                'peak':  int(peak),
            })

    return episodes

# ── Match episodes ────────────────────────────────────────────────────────────

def match_episodes(ep_ref, ep_cand, max_shift=60):
    """
    Match retrograde episodes between reference and candidate series.
    Each reference episode is matched to the nearest candidate episode
    within max_shift days.
    Returns list of (ref_ep, cand_ep, shift_start, shift_peak, shift_end).
    """
    matches = []
    used    = set()

    for r in ep_ref:
        best_match = None
        best_dist  = max_shift + 1

        for ci, c in enumerate(ep_cand):
            if ci in used:
                continue
            dist = abs(r['start'] - c['start'])
            if dist < best_dist:
                best_dist  = dist
                best_match = (ci, c)

        if best_match is not None:
            ci, c = best_match
            used.add(ci)
            matches.append({
                'ref_start':  r['start'],
                'ref_peak':   r['peak'],
                'ref_end':    r['end'],
                'ref_dur':    r['dur'],
                'cand_start': c['start'],
                'cand_peak':  c['peak'],
                'cand_end':   c['end'],
                'cand_dur':   c['dur'],
                'shift_start': r['start'] - c['start'],  # days
                'shift_peak':  r['peak']  - c['peak'],
                'shift_end':   r['end']   - c['end'],
                'dur_diff':    r['dur']   - c['dur'],
            })

    return matches

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  Retrograde synchronisation analysis")
    print(f"  T* = {T_STAR} days ({T_STAR/365.25:.4f} years)")
    print("=" * 60)

    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START
    print(f"\nLoading ephemeris...")
    eph = load(EPHEMERIS)
    print("Ephemeris loaded.")

    t_ref = ts.utc(REF_YEAR, REF_MONTH, REF_DAY, 12)
    jd_ref = int(t_ref.tdb)
    print(f"Reference epoch: {REF_YEAR}/{REF_MONTH}/{REF_DAY}  (JD {jd_ref})")

    # Load/build geocentric cache
    positions = load_or_build_geo(ts, eph, jd_ref, SERIES_DAYS, T_STAR)

    all_episode_rows = []
    stats_rows       = []

    print()
    for name in PLANET_NAMES:
        lon_ref, lon_cand = positions[name]
        min_days = MIN_RETRO_DAYS[name]

        ep_ref  = find_retrogrades(lon_ref,  min_days)
        ep_cand = find_retrogrades(lon_cand, min_days)

        print(f"\n{'─'*60}")
        print(f"  {name}  (min duration: {min_days} days)")
        print(f"  Retrograde episodes — reference series: {len(ep_ref)}")
        print(f"  Retrograde episodes — candidate series: {len(ep_cand)}")

        max_shift = 400 if name == 'Neptune' else 60
        matches = match_episodes(ep_ref, ep_cand, max_shift=max_shift)
        
        print(f"  Matched pairs: {len(matches)}")

        if not matches:
            continue

        shifts_start = np.array([m['shift_start'] for m in matches])
        shifts_peak  = np.array([m['shift_peak']  for m in matches])
        shifts_end   = np.array([m['shift_end']   for m in matches])
        dur_diffs    = np.array([m['dur_diff']     for m in matches])

        print(f"\n  Shift statistics (reference − candidate, days):")
        print(f"  {'':10} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
        print(f"  {'Start':10} {shifts_start.mean():>8.2f} "
              f"{shifts_start.std():>8.2f} "
              f"{shifts_start.min():>8.1f} {shifts_start.max():>8.1f}")
        print(f"  {'Peak':10} {shifts_peak.mean():>8.2f} "
              f"{shifts_peak.std():>8.2f} "
              f"{shifts_peak.min():>8.1f} {shifts_peak.max():>8.1f}")
        print(f"  {'End':10} {shifts_end.mean():>8.2f} "
              f"{shifts_end.std():>8.2f} "
              f"{shifts_end.min():>8.1f} {shifts_end.max():>8.1f}")
        print(f"  {'Duration':10} {dur_diffs.mean():>8.2f} "
              f"{dur_diffs.std():>8.2f} "
              f"{dur_diffs.min():>8.1f} {dur_diffs.max():>8.1f}")

        # Convert to hours for planets with small shifts
        print(f"\n  Peak shift in hours: "
              f"mean={shifts_peak.mean()*24:.1f}h  "
              f"std={shifts_peak.std()*24:.1f}h  "
              f"max={abs(shifts_peak).max()*24:.1f}h")

        # Save episode rows
        for m in matches:
            row = {'planet': name}
            row.update(m)
            all_episode_rows.append(row)

        stats_rows.append({
            'planet':          name,
            'n_ref':           len(ep_ref),
            'n_cand':          len(ep_cand),
            'n_matched':       len(matches),
            'mean_shift_start_d': round(shifts_start.mean(), 3),
            'std_shift_start_d':  round(shifts_start.std(),  3),
            'mean_shift_peak_d':  round(shifts_peak.mean(),  3),
            'std_shift_peak_d':   round(shifts_peak.std(),   3),
            'mean_shift_end_d':   round(shifts_end.mean(),   3),
            'std_shift_end_d':    round(shifts_end.std(),    3),
            'mean_shift_peak_h':  round(shifts_peak.mean()*24, 1),
            'std_shift_peak_h':   round(shifts_peak.std()*24,  1),
            'max_abs_shift_peak_h': round(abs(shifts_peak).max()*24, 1),
            'mean_dur_diff_d':    round(dur_diffs.mean(), 3),
        })

    # Save CSV files
    if all_episode_rows:
        fields = list(all_episode_rows[0].keys())
        with open('retrograde_sync.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(all_episode_rows)
        print(f"\n  -> retrograde_sync.csv")

    if stats_rows:
        fields = list(stats_rows[0].keys())
        with open('retrograde_stats.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(stats_rows)
        print(f"  -> retrograde_stats.csv")

    # Summary table
    print(f"\n{'═'*60}")
    print(f"  Summary: retrograde peak shift at T* = {T_STAR} days")
    print(f"{'═'*60}")
    print(f"  {'Planet':<10} {'Episodes':>9} {'Mean shift':>12} "
          f"{'Std shift':>11} {'Max shift':>10}")
    print(f"  {'':10} {'(matched)':>9} {'(hours)':>12} "
          f"{'(hours)':>11} {'(hours)':>10}")
    print(f"  {'─'*55}")
    for r in stats_rows:
        print(f"  {r['planet']:<10} {r['n_matched']:>9} "
              f"{r['mean_shift_peak_h']:>+12.1f} "
              f"{r['std_shift_peak_h']:>11.1f} "
              f"{r['max_abs_shift_peak_h']:>10.1f}")
    print(f"{'═'*60}")


if __name__ == '__main__':
    main()
