#!/usr/bin/env python3
"""
helio1151.py
============
Heliocentric analysis of the multi-planet quasi-commensurability.

Question:
  What is the time interval T that minimises the mean simultaneous
  angular displacement of the five naked-eye planets from their
  heliocentric positions T days earlier?

Bodies: Mercury, Venus, Mars, Jupiter, Saturn (5 naked-eye planets).
Frame:  Heliocentric ecliptic longitude, J2000.
Source: DE441-part1 ephemeris via Skyfield.

Analyses performed:
  1. Full cycle search: exhaustive search over [-4000, +1500] years
  2. Statistical significance of the global minimum
  3. Temporal stability across 12 reference epochs
  4. Per-planet breakdown at T* 
  5. Theoretical sidereal residues at T*
  6. Secondary minima characterisation
  7. Daily offset extraction for figures (fast + slow planet panels)
  8. Polar snapshot data for Figure C (5 epochs)

Output files:
  helio_cache.pkl          position cache (heliocentric)
  helio_results.csv        full search results
  helio_stats.csv          per-planet statistics at T*
  helio_panel_fast.csv     Mercury, Venus, Mars offsets — 5yr daily
  helio_panel_slow.csv     Jupiter, Saturn offsets — 100yr weekly
  helio_polar.csv          polar snapshots at 5 epochs

Ephemeris:
  de441-part1.bsp
  Download: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp
  Covers: ~-13200 to ~+1550 CE

Usage:
  pip install skyfield numpy scipy
  python helio1151.py

Author: Carlos Baiget Orts — asinfreedom@gmail.com
"""

import os
import pickle
import csv
import lzma
import numpy as np
from skyfield.api import load, GREGORIAN_START
from skyfield.framelib import ecliptic_J2000_frame

# ── Configuration ─────────────────────────────────────────────────────────────

# de441_part-1.bsp: covers 13201 BCE to 1969 CE (~1.5 GB)
# With ref year 0 and search ±1300yr, coverage needed is -1300 to +1300 CE
# — well within part-1. No need for part-2.
# Skyfield downloads automatically if not present locally.
EPHEMERIS    = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp'
CACHE_FILE   = 'helio6_de441_cache.pkl.xz'  # 6 pre-Uranian planets, de441-part1, LZMA compressed

# Reference epoch (arbitrary — stability test confirms independence)
REF_YEAR, REF_MONTH, REF_DAY = 0, 6, 15  # year 0: search ±1300yr spans -1300 to +1300 CE

# Series length: 100 Julian years
SERIES_DAYS  = 36525

# Search range (years relative to reference epoch)
# de441-part1 covers up to ~1550 CE; ref year 1 + 1500 = 1501 CE — safe
SEARCH_MIN   = -1300  # symmetric around year 0: spans -1300 to +1300 CE
SEARCH_MAX   = +1300  # both within de441-part1 coverage (up to 1969 CE)
STEP_YEARS   = 1

# Bodies: 6 planets visible to the naked eye from outside the solar system
# Earth is included — from a heliocentric perspective it is a planet like the others
PLANET_NAMES    = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn']
PLANET_TARGETS  = [
    'mercury barycenter',
    'venus barycenter',
    'earth',               # Earth barycenter ≈ Earth (Moon mass < 0.01%)
    'mars barycenter',
    'jupiter barycenter',
    'saturn barycenter',
]
N_PLANETS = len(PLANET_NAMES)

# Known optimal value (used for breakdown and figure extraction)
T_STAR = 420403   # days — emerges from search, set here for convenience

# ── Ephemeris and timescale ───────────────────────────────────────────────────

def init_timescale():
    """Initialise Skyfield timescale only — no ephemeris download."""
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START
    return ts

def load_ephemeris():
    """Load (and download if necessary) the DE441 ephemeris."""
    print("Loading ephemeris (downloads if not present locally)...")
    eph = load(EPHEMERIS)
    print("Ephemeris loaded.")
    return eph

# ── Cache: heliocentric ecliptic longitudes ───────────────────────────────────

def build_cache(ts, eph, jd_start, jd_end):
    """
    Compute daily heliocentric ecliptic longitudes (J2000) for the
    5 naked-eye planets over [jd_start, jd_end].

    Observer: Sun (heliocentric).
    Frame: ecliptic J2000.

    Returns: (jd_base:int, data:np.ndarray shape (N_days, 5) float32)
    """
    print(f"\nBuilding heliocentric position cache...")
    print(f"  Range: JD {int(jd_start)} to {int(jd_end)}")
    print(f"  Bodies: {', '.join(PLANET_NAMES)}")
    print(f"  (This may take several minutes — saved to {CACHE_FILE})\n")

    sun     = eph['sun']
    targets = [eph[t] for t in PLANET_TARGETS]
    n_days  = int(jd_end - jd_start) + 1
    jd_base = int(jd_start)
    data    = np.zeros((n_days, N_PLANETS), dtype=np.float32)

    batch = 500
    for start in range(0, n_days, batch):
        end  = min(start + batch, n_days)
        jds  = np.arange(jd_base + start, jd_base + end, dtype=float)
        t_b  = ts.tt_jd(jds)

        for pi, tgt in enumerate(targets):
            # Heliocentric position: observe planet from the Sun
            astrometric = sun.at(t_b).observe(tgt)
            lat, lon, _ = astrometric.frame_latlon(ecliptic_J2000_frame)
            data[start:end, pi] = lon.degrees.astype(np.float32)

        pct = round(end / n_days * 100, 1)
        print(f"  {pct}% ({end}/{n_days} days)", end='\r')

    print(f"\n  Cache built: {n_days} days × {N_PLANETS} planets.")
    return jd_base, data


def load_or_build_cache(ts, jd_ref, series_days, search_min_y, search_max_y, eph=None):
    margin  = series_days + 10
    jd_start = jd_ref + search_min_y * 365.25 - margin
    jd_end   = jd_ref + search_max_y * 365.25 + margin

    if os.path.exists(CACHE_FILE):
        print(f"Loading cache from {CACHE_FILE}...")
        with lzma.open(CACHE_FILE, 'rb') as f:
            saved = pickle.load(f)
        ok = (saved['jd_base'] <= jd_start and
              saved['jd_base'] + len(saved['data']) >= jd_end)
        if ok:
            print("Cache OK.")
            return saved['jd_base'], saved['data']
        print("Cache range insufficient — rebuilding.")

    if eph is None:
        eph = load_ephemeris()
    jd_base, data = build_cache(ts, eph, jd_start, jd_end)
    with lzma.open(CACHE_FILE, 'wb') as f:
        pickle.dump({'jd_base': jd_base, 'data': data}, f)
    print(f"Cache saved to {CACHE_FILE}\n")
    return jd_base, data

# ── Core metric ───────────────────────────────────────────────────────────────

def get_series(data, jd_base, jd_start, length):
    idx = int(jd_start) - jd_base
    return data[idx:idx + length]   # shape (length, N_PLANETS)


def angular_dist(a, b):
    """Unsigned circular angular distance, degrees, result in [0, 180]."""
    d = np.abs(a.astype(np.float64) - b.astype(np.float64)) % 360.0
    return np.where(d > 180.0, 360.0 - d, d)


def signed_diff(a, b):
    """Signed circular difference a − b, degrees, result in (−180, +180]."""
    d = (a.astype(np.float64) - b.astype(np.float64)) % 360.0
    return np.where(d > 180.0, d - 360.0, d)


def score(series_ref, series_cand):
    """
    S(T) = mean(daily_means) + std(daily_means)
    where daily_mean_i = mean angular distance across all planets on day i.
    """
    diff        = angular_dist(series_ref, series_cand)   # (N, 5)
    daily_means = diff.mean(axis=1)                        # (N,)
    return float(daily_means.mean()), float(daily_means.std())

# ── Analysis 1: full cycle search ────────────────────────────────────────────

def run_search(jd_ref, jd_base, data, series_days, search_min, search_max, step):
    s_ref   = get_series(data, jd_base, jd_ref, series_days)
    results = []
    cands   = list(range(search_min, search_max + 1, step))
    total   = len(cands)

    for i, dy in enumerate(cands):
        if i % 200 == 0:
            print(f"  Searching... {round(i/total*100,1)}%", end='\r')
        jd_c = jd_ref + dy * 365.25
        try:
            s_c = get_series(data, jd_base, jd_c, series_days)
            if len(s_c) < series_days:
                continue
            m, s = score(s_ref, s_c)
            results.append({
                'delta_years': dy,
                'delta_days':  round(dy * 365.25),
                'mean_deg':    round(m, 4),
                'std_deg':     round(s, 4),
                'score':       round(m + s, 4),
            })
        except Exception:
            continue

    results.sort(key=lambda r: r['score'])
    print(f"\n  Search complete. {len(results)} candidates evaluated.")
    return results

# ── Analysis 2: statistical significance ─────────────────────────────────────

def statistical_significance(results, best_score):
    scores = np.array([r['score'] for r in results if r['delta_years'] != 0])
    mean_s = scores.mean()
    std_s  = scores.std()
    z      = (mean_s - best_score) / std_s
    pct    = (scores > best_score).mean() * 100

    print(f"\n{'─'*55}")
    print(f"Statistical significance")
    print(f"{'─'*55}")
    print(f"  Candidates:          {len(scores)}")
    print(f"  Mean score (all):    {mean_s:.3f}°")
    print(f"  Std  score (all):    {std_s:.3f}°")
    print(f"  Median score:        {np.median(scores):.3f}°")
    print(f"  Best score (T*):     {best_score:.3f}°")
    print(f"  Max score:           {scores.max():.3f}°")
    print(f"  Z-score:             {z:.2f}σ")
    print(f"  Percentile:          {pct:.1f}% of candidates are worse")

    bins = [0,15,18,20,22,24,26,28,30,35,40,50,60,80,110]
    counts, _ = np.histogram(scores, bins=bins)
    print(f"\n  Score distribution:")
    for i in range(len(counts)):
        bar = '█' * int(counts[i] / len(scores) * 150)
        print(f"    {bins[i]:>4}–{bins[i+1]:<4}°: {counts[i]:>5}  {bar}")
    print()
    return z, pct

# ── Analysis 3: temporal stability ───────────────────────────────────────────

def stability_test(ts, jd_base, data, series_days, t_star, n=12):
    print(f"\n{'─'*55}")
    print(f"Temporal stability (T* = {t_star} days, {n} reference epochs)")
    print(f"{'─'*55}")
    sc_list = []
    # Valid epochs: need epoch - T* >= cache_start and epoch + series <= cache_end
    # With cache ~-1300 to +1300 and T*=1151yr: epoch >= -1300+1151 = -149
    # Use epochs from -100 to +1100 to stay safely within cache
    for i in range(n):
        year = -100 + i * 110
        t    = ts.utc(year, 6, 15, 12)
        jd   = int(t.tdb)
        try:
            sr = get_series(data, jd_base, jd,          series_days)
            sc = get_series(data, jd_base, jd - t_star, series_days)
            if len(sr) < series_days or len(sc) < series_days:
                continue
            m, s = score(sr, sc)
            sc_list.append(m + s)
            print(f"  Ref year {year:>6}: mean={m:.3f}°  std={s:.3f}°  "
                  f"score={m+s:.3f}°")
        except Exception as e:
            print(f"  Ref year {year:>6}: skipped ({e})")
    if sc_list:
        print(f"\n  Score range:  {min(sc_list):.3f}° – {max(sc_list):.3f}°")
        print(f"  Mean score:   {np.mean(sc_list):.3f}°")
        print(f"  Std of scores:{np.std(sc_list):.4f}°  ← lower = more stable")
    print()

# ── Analysis 4: per-planet breakdown ─────────────────────────────────────────

def per_planet_breakdown(jd_ref, jd_base, data, series_days, t_star):
    sr = get_series(data, jd_base, jd_ref,          series_days)
    sc = get_series(data, jd_base, jd_ref - t_star, series_days)

    diff   = angular_dist(sr, sc)    # (N, 5) unsigned
    sdiff  = signed_diff(sr, sc)     # (N, 5) signed

    mean_abs  = diff.mean(axis=0)
    std_abs   = diff.std(axis=0)
    mean_sign = sdiff.mean(axis=0)

    print(f"\n{'─'*72}")
    print(f"Per-planet breakdown at T* = {t_star} days "
          f"({t_star/365.25:.4f} years)")
    print(f"{'─'*72}")
    print(f"{'Planet':<10} {'Mean|Δ|(°)':>12} {'Std|Δ|(°)':>11} "
          f"{'Signed mean(°)':>16} {'Range(°)':>10}")
    print(f"{'─'*72}")

    stats = []
    for k, name in enumerate(PLANET_NAMES):
        rng = float(sdiff[:, k].max() - sdiff[:, k].min())
        print(f"{name:<10} {mean_abs[k]:>12.3f} {std_abs[k]:>11.3f} "
              f"{mean_sign[k]:>+16.3f} {rng:>10.3f}")
        stats.append({
            'planet':      name,
            'mean_abs':    round(float(mean_abs[k]),  4),
            'std_abs':     round(float(std_abs[k]),   4),
            'mean_signed': round(float(mean_sign[k]), 4),
            'min_signed':  round(float(sdiff[:,k].min()), 4),
            'max_signed':  round(float(sdiff[:,k].max()), 4),
            'range':       round(rng, 4),
        })

    total_mean = mean_abs.mean()
    total_std  = mean_abs.mean(axis=0)
    dm = diff.mean(axis=1)
    print(f"{'─'*72}")
    print(f"{'TOTAL':<10} {total_mean:>12.3f}  "
          f"score = {dm.mean():.3f}° + {dm.std():.3f}° = "
          f"{dm.mean()+dm.std():.3f}°")
    print()

    with open('helio_stats.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(stats[0].keys()))
        writer.writeheader()
        writer.writerows(stats)
    print("  -> helio_stats.csv")
    return stats

# ── Analysis 5: theoretical sidereal residues ─────────────────────────────────

def sidereal_residues(t_star):
    """
    For each planet, compute how many degrees short of an integer number
    of complete orbits T* days represents.
    """
    sidereal = {
        'Mercury':   87.9691,
        'Venus':    224.7008,
        'Earth':    365.2500,
        'Mars':     686.9710,
        'Jupiter': 4332.589,
        'Saturn':  10759.22,
    }
    print(f"\n{'─'*68}")
    print(f"Theoretical sidereal residues at T* = {t_star} days "
          f"({t_star/365.25:.4f} yr)")
    print(f"{'─'*68}")
    print(f"{'Planet':<10} {'Sid. period(d)':>15} {'Orbits in T*':>14} "
          f"{'Frac. res.':>12} {'Residue(°)':>12}")
    print(f"{'─'*68}")

    residues = []
    for name in PLANET_NAMES:
        p     = sidereal[name]
        cyc   = t_star / p
        frac  = cyc - int(cyc)
        if frac > 0.5:
            frac -= 1
        deg   = frac * 360
        residues.append(deg)
        print(f"{name:<10} {p:>15.4f} {cyc:>14.4f} "
              f"{frac:>+12.6f} {deg:>+12.3f}°")

    print(f"{'─'*68}")
    print(f"  Mean |residue|: {np.mean(np.abs(residues)):.3f}°")
    print(f"  Note: Earth year residue = "
          f"{((t_star/365.25 - int(t_star/365.25))*360):.3f}°  "
          f"(T* ≈ integer number of years → helio ≈ geo)")
    print()
    return residues

# ── Analysis 6: secondary minima ─────────────────────────────────────────────

def secondary_minima(results, top_n=15):
    pjs = 19.859   # Jupiter-Saturn conjunction period (years)
    print(f"\n{'─'*70}")
    print(f"Top {top_n} best cycles")
    print(f"{'─'*70}")
    print(f"{'Rank':>4} {'ΔYr':>8} {'ΔDays':>9} "
          f"{'Mean°':>8} {'Std°':>7} {'Score°':>8}  Relation")
    print(f"{'─'*70}")
    shown = 0
    for r in results:
        if r['delta_years'] == 0:
            continue
        dy  = r['delta_years']
        rat = dy / 1151 if dy != 0 else 0
        # Check proximity to simple fractions of 1151
        fracs = [(n,d) for d in range(1,5) for n in range(-5*d, 6*d) if n!=0]
        best  = min(fracs, key=lambda f: abs(rat - f[0]/f[1]))
        brat  = best[0]/best[1]
        if abs(rat - brat) < 0.04:
            rel = f"{best[0]}/{best[1]} × 1151"
        else:
            # Check JS conjunction multiples
            js_mult = round(abs(dy) / pjs)
            if abs(abs(dy) - js_mult * pjs) < 8:
                rel = f"≈ {js_mult} × P_JS"
            else:
                rel = f"~{rat:.2f} × 1151"
        print(f"{shown+1:>4} {dy:>8} {r['delta_days']:>9} "
              f"{r['mean_deg']:>8.3f} {r['std_deg']:>7.3f} "
              f"{r['score']:>8.3f}  {rel}")
        shown += 1
        if shown >= top_n:
            break
    print()

# ── Analysis 7: daily offsets for figures ────────────────────────────────────

def extract_figure_data(jd_ref, jd_base, data, series_days, t_star):
    sr = get_series(data, jd_base, jd_ref,          series_days)
    sc = get_series(data, jd_base, jd_ref - t_star, series_days)
    sd = signed_diff(sr, sc)   # shape (series_days, 5)
    years = np.arange(series_days) / 365.25

    # Fast panel: Mercury(0), Venus(1), Mars(2) — 5 years, daily
    fast_days = int(5 * 365.25)
    with open('helio_panel_fast.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['day', 'year', 'Mercury', 'Venus', 'Earth', 'Mars'])
        for i in range(fast_days):
            w.writerow([i, round(float(years[i]),5),
                        round(float(sd[i,0]),3),
                        round(float(sd[i,1]),3),
                        round(float(sd[i,2]),3),
                        round(float(sd[i,3]),3)])
    print("  -> helio_panel_fast.csv  (Mercury, Venus, Earth, Mars — 5yr daily)")

    # Slow panel: Jupiter(3), Saturn(4) — 100 years, weekly
    with open('helio_panel_slow.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['day', 'year', 'Jupiter', 'Saturn'])
        for i in range(0, series_days, 7):
            w.writerow([i, round(float(years[i]),4),
                        round(float(sd[i,4]),3),
                        round(float(sd[i,5]),3)])
    print("  -> helio_panel_slow.csv  (Jupiter, Saturn — 100yr weekly)")

# ── Analysis 8: polar snapshot data ──────────────────────────────────────────

def extract_polar_data(ts, jd_base, data, t_star):
    snapshot_years = [50, 200, 400, 600, 800]  # epoch-1151 range: -1101 to -351 CE — safely within cache
    rows = []
    for yr in snapshot_years:
        t   = ts.utc(yr, 6, 15, 12)
        jd  = int(t.tdb)
        try:
            pr = get_series(data, jd_base, jd,          1)[0]
            pc = get_series(data, jd_base, jd - t_star, 1)[0]
            for k, name in enumerate(PLANET_NAMES):
                off = signed_diff(np.array([pr[k]]),
                                  np.array([pc[k]]))[0]
                rows.append({
                    'epoch': yr,
                    'planet': name,
                    'lon_ref':  round(float(pr[k]), 3),
                    'lon_cand': round(float(pc[k]), 3),
                    'offset':   round(float(off),   3),
                })
        except Exception as e:
            print(f"  Polar snapshot year {yr}: skipped ({e})")

    with open('helio_polar.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['epoch','planet',
                                          'lon_ref','lon_cand','offset'])
        w.writeheader()
        w.writerows(rows)
    print("  -> helio_polar.csv  (5 epochs × 6 planets)")

    # Print offset summary per epoch
    print(f"\n  Offset summary (signed, degrees):")
    print(f"  {'Planet':<10}", end='')
    for yr in snapshot_years:
        print(f"  {yr:>6}CE", end='')
    print()
    for name in PLANET_NAMES:
        print(f"  {name:<10}", end='')
        for yr in snapshot_years:
            val = next(r['offset'] for r in rows
                       if r['epoch']==yr and r['planet']==name)
            print(f"  {val:>+8.2f}°", end='')
        print()
    print()


def series_length_convergence(jd_ref, jd_base, data, t_star):
    """
    Demonstrate that the result is stable across series lengths.
    For each series length in LENGTHS:
      - Run the full search and record the global minimum and its score
      - Record the score of T* specifically
      - Record the rank of T* among all candidates
    Shows that T* emerges as the global minimum even with short series,
    and that scores converge quickly as series length increases.
    """
    # Series lengths to test: from 1 year to 100 years
    LENGTHS = [
        (1,    365),
        (2,    730),
        (5,   1826),
        (10,  3652),
        (20,  7305),
        (50, 18262),
        (100, 36525),
    ]

    print(f"\n{'─'*72}")
    print(f"Series length convergence analysis")
    print(f"T* = {t_star} days — does it emerge as global minimum for short series?")
    print(f"{'─'*72}")
    print(f"{'Length':>8} {'N days':>8} {'Best T (yr)':>12} {'Best score':>12} "
          f"{'T* score':>10} {'T* rank':>8} {'2nd best':>10}")
    print(f"{'─'*72}")

    conv_rows = []
    s_ref_full = get_series(data, jd_base, jd_ref, 36525)

    for label, n_days in LENGTHS:
        s_ref = s_ref_full[:n_days]
        results_l = []

        for dy in range(-1300, 1301, 1):
            if dy == 0:
                continue
            jd_c = jd_ref + dy * 365.25
            try:
                s_c = get_series(data, jd_base, jd_c, n_days)
                if len(s_c) < n_days:
                    continue
                m, s = score(s_ref, s_c)
                results_l.append((dy, round(m+s, 4)))
            except Exception:
                continue

        results_l.sort(key=lambda x: x[1])
        best_dy, best_sc = results_l[0]
        second_sc = results_l[1][1] if len(results_l) > 1 else float('nan')

        # Score and rank of T*
        tstar_scores = [(dy, sc) for dy, sc in results_l if abs(abs(dy) - 1151) <= 1]
        if tstar_scores:
            tstar_sc = min(s for _, s in tstar_scores)
            tstar_rank = next(i+1 for i, (dy, sc) in enumerate(results_l)
                              if abs(abs(dy) - 1151) <= 1)
        else:
            tstar_sc = float('nan')
            tstar_rank = -1

        print(f"{label:>7}yr {n_days:>8} {best_dy:>+12} {best_sc:>12.3f} "
              f"{tstar_sc:>10.3f} {tstar_rank:>8} {second_sc:>10.3f}")

        conv_rows.append({
            'series_years':  label,
            'series_days':   n_days,
            'best_dy':       best_dy,
            'best_score':    best_sc,
            'tstar_score':   tstar_sc,
            'tstar_rank':    tstar_rank,
            'second_score':  second_sc,
        })

    print(f"{'─'*72}")
    print()

    # Save to CSV for figure
    import csv as csv_mod
    with open('helio_convergence.csv', 'w', newline='') as f:
        w = csv_mod.DictWriter(f, fieldnames=list(conv_rows[0].keys()))
        w.writeheader()
        w.writerows(conv_rows)
    print("  -> helio_convergence.csv")
    return conv_rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("  Heliocentric Quasi-Commensurability Analysis")
    print("  Pre-Uranian planets: Mercury · Venus · Earth · Mars · Jupiter · Saturn")
    print("=" * 65)
    print()

    ts = init_timescale()

    t_ref = ts.utc(REF_YEAR, REF_MONTH, REF_DAY, 12)
    jd_ref = int(t_ref.tdb)

    print(f"\nReference epoch: {REF_YEAR}/{REF_MONTH}/{REF_DAY}  (JD {jd_ref})")
    print(f"Series length:   {SERIES_DAYS} days ({SERIES_DAYS/365.25:.1f} yr)")
    print(f"Search range:    {SEARCH_MIN} to {SEARCH_MAX} years")
    print(f"Candidates:      {SEARCH_MAX - SEARCH_MIN + 1}")
    print()

    # Load cache — ephemeris downloaded only if cache is absent or insufficient
    jd_base, data = load_or_build_cache(
        ts, jd_ref, SERIES_DAYS, SEARCH_MIN, SEARCH_MAX
    )

    # ── 1. Full search ────────────────────────────────────────────────────────
    print("\n[1/9] Full cycle search...")
    results = run_search(
        jd_ref, jd_base, data,
        SERIES_DAYS, SEARCH_MIN, SEARCH_MAX, STEP_YEARS
    )
    best = next(r for r in results if r['delta_years'] != 0)
    print(f"\n{'═'*50}")
    print(f"  GLOBAL MINIMUM:")
    print(f"  T*      = {best['delta_days']:,} days")
    print(f"  T* / yr = {best['delta_days']/365.25:.6f} years")
    print(f"  Mean δ  = {best['mean_deg']:.4f}°")
    print(f"  Std  δ  = {best['std_deg']:.4f}°")
    print(f"  Score   = {best['score']:.4f}°")
    print(f"{'═'*50}")

    with open('helio_results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f,
            fieldnames=['delta_years','delta_days','mean_deg','std_deg','score'])
        w.writeheader()
        w.writerows(results)
    print("  -> helio_results.csv")

    # ── 2. Statistical significance ───────────────────────────────────────────
    print("\n[2/9] Statistical significance...")
    statistical_significance(results, best['score'])

    # ── 3. Temporal stability ─────────────────────────────────────────────────
    print("\n[3/9] Temporal stability...")
    stability_test(ts, jd_base, data, SERIES_DAYS, T_STAR)

    # ── 4. Per-planet breakdown ───────────────────────────────────────────────
    print("\n[4/9] Per-planet breakdown at T* = 420,403 days...")
    per_planet_breakdown(jd_ref, jd_base, data, SERIES_DAYS, T_STAR)

    # ── 5. Theoretical sidereal residues ──────────────────────────────────────
    print("\n[5/9] Theoretical sidereal residues...")
    sidereal_residues(T_STAR)

    # ── 6. Secondary minima ───────────────────────────────────────────────────
    print("\n[6/9] Secondary minima...")
    secondary_minima(results, top_n=15)

    # ── 7. Figure data: daily offsets ─────────────────────────────────────────
    print("\n[7/9] Extracting figure data (panel plots)...")
    extract_figure_data(jd_ref, jd_base, data, SERIES_DAYS, T_STAR)

    # ── 8. Figure data: polar snapshots ───────────────────────────────────────
    print("\n[8/9] Extracting polar snapshot data...")
    extract_polar_data(ts, jd_base, data, T_STAR)

    # ── 9. Series length convergence ──────────────────────────────────────────
    print("\n[9/9] Series length convergence analysis...")
    series_length_convergence(jd_ref, jd_base, data, T_STAR)

    print("\n" + "=" * 65)
    print("  All analyses complete. Output files:")
    print("    helio6_de441_cache.pkl.xz  position cache (LZMA compressed)")
    print("    helio_results.csv        full search — scatter plot data")
    print("    helio_stats.csv          per-planet statistics at T*")
    print("    helio_panel_fast.csv     Mercury, Venus, Earth, Mars — Fig A top")
    print("    helio_panel_slow.csv     Jupiter, Saturn — Fig A bottom")
    print("    helio_polar.csv          5 epochs x 6 planets — Fig C")
    print("    helio_convergence.csv    series length convergence — Fig D")
    print("=" * 65)


if __name__ == '__main__':
    main()
