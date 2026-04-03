#!/usr/bin/env python3
"""
make_figures.py
===============
Generates publication-quality figures for the 1151-year heliocentric
quasi-commensurability paper.

Figures produced:
  fig1_scatter.pdf    — Score S(T) vs candidate interval (all 5501 candidates)
  fig2_histogram.pdf  — Score distribution histogram
  fig3_panels.pdf     — Angular offset per planet over time (2 panels)
  fig4_polar.pdf      — Polar snapshots at 5 epochs with per-planet radii

Input files (from helio1151.py):
  helio_results.csv
  helio_panel_fast.csv
  helio_panel_slow.csv
  helio_polar.csv

Usage:
  pip install matplotlib numpy
  python make_figures.py

Author: Carlos Baiget Orts -- asinfreedom@gmail.com
"""

import numpy as np
import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family':       'serif',
    'font.size':         10,
    'axes.titlesize':    11,
    'axes.labelsize':    10,
    'xtick.labelsize':    9,
    'ytick.labelsize':    9,
    'legend.fontsize':    9,
    'figure.dpi':        150,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         True,
    'grid.alpha':        0.3,
    'grid.linewidth':    0.5,
    'lines.linewidth':   1.2,
})

COLORS = {
    'Mercury': '#888888',
    'Venus':   '#E07B39',
    'Earth':   '#4A90D9',
    'Mars':    '#C0392B',
    'Jupiter': '#8E44AD',
    'Saturn':  '#16A085',
}

POLAR_RADIUS = {
    'Mercury': 0.35,
    'Venus':   0.50,
    'Earth':   0.62,
    'Mars':    0.72,
    'Jupiter': 0.84,
    'Saturn':  0.95,
}


def load_csv(path):
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) if k not in ('planet',) else v
                         for k, v in row.items()})
    return rows


def load_results(path='helio_results.csv'):
    rows = load_csv(path)
    dy = np.array([r['delta_years'] for r in rows])
    sc = np.array([r['score']       for r in rows])
    return dy, sc


def load_polar(path='helio_polar.csv'):
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                'epoch':    int(float(row['epoch'])),
                'planet':   row['planet'],
                'lon_ref':  float(row['lon_ref']),
                'lon_cand': float(row['lon_cand']),
                'offset':   float(row['offset']),
            })
    return rows


def fig1_scatter():
    dy, sc = load_results()
    fig, ax = plt.subplots(figsize=(8, 4))

    mask_zero = (dy == 0)
    mask_best = (np.abs(dy) == 1151)
    mask_rest = ~mask_zero & ~mask_best

    ax.scatter(dy[mask_rest], sc[mask_rest],
               s=1.5, c='#1a3a5c', alpha=1, linewidths=0,
               label='Candidates')
    # Draw the actual data point (same style as others but slightly larger)
    ax.scatter(dy[mask_best], sc[mask_best],
               s=8, c='#1a3a5c', alpha=1, linewidths=0, zorder=5)
    # Draw red circle around it
    ax.scatter(dy[mask_best], sc[mask_best],
               s=40, facecolors='none', edgecolors='#C0392B',
               linewidths=1.8, zorder=6,
               label=r'$T^* = \pm\,1151\,\mathrm{yr}$')

    # Annotations: -1151 label to the right, +1151 label to the left
    for sign, lbl, xoff, yoff, ha in [
        (-1, r'$-1151\,\mathrm{yr}$',  200,  18, 'left'),
        ( 1, r'$+1151\,\mathrm{yr}$', -200,  18, 'right'),
    ]:
        idx = np.where(dy == sign * 1151)[0]
        if len(idx) == 0:
            continue
        i = idx[0]
        ax.annotate(lbl,
                    xy=(dy[i], sc[i]),
                    xytext=(dy[i] + xoff, sc[i] + yoff),
                    arrowprops=dict(arrowstyle='->', color='#C0392B', lw=1.2,
                                   shrinkA=8, shrinkB=4),
                    fontsize=9, color='#C0392B', ha=ha,
                    bbox=dict(fc='white', ec='none', pad=1))

    ax.set_xlabel(r'$\Delta t$ (years)')
    ax.set_ylabel(r'Score $S(\Delta t)$ (degrees)')
    ax.set_xlim(-1400, 1400)
    ax.set_ylim(0, 150)
    ax.legend(loc='upper right', framealpha=0.9, markerscale=2)
    fig.tight_layout()
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig1_scatter.pdf', bbox_inches='tight')
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig1_scatter.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  -> fig1_scatter.pdf / .png")


def fig2_histogram():
    dy, sc = load_results()
    sc_nz = sc[dy != 0]
    fig, ax = plt.subplots(figsize=(7, 3.8))

    bins = np.arange(0, 145, 3)
    n, edges, patches = ax.hist(sc_nz, bins=bins,
                                color='#4A90D9', alpha=0.75,
                                edgecolor='white', linewidth=0.3)
    for patch, left in zip(patches, edges[:-1]):
        if left < 18:
            patch.set_facecolor('#C0392B')
            patch.set_alpha(1.0)

    best = sc_nz.min()
    ax.axvline(best, color='#C0392B', lw=1.2, linestyle='--')
    ax.text(best + 1.5, n.max() * 0.85,
            r'$T^*$' + f' = {best:.1f}' + r'$\degree$',
            color='#C0392B', fontsize=8.5)

    z = (sc_nz.mean() - best) / sc_nz.std()
    stats_text = (r'$n$ = ' + f'{len(sc_nz):,}' + '\n'
                  + r'Mean = ' + f'{sc_nz.mean():.1f}' + r'$\degree$' + '\n'
                  + r'$z$ = ' + f'{z:.2f}' + r'$\,\sigma$')
    ax.text(0.97, 0.95, stats_text,
            transform=ax.transAxes, ha='right', va='top',
            fontsize=8.5, bbox=dict(fc='white', ec='#cccccc', lw=0.5, pad=4))

    ax.set_xlabel(r'Score $S(\Delta t)$ (degrees)')
    ax.set_ylabel('Number of candidates')
    fig.tight_layout()
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig2_histogram.pdf', bbox_inches='tight')
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig2_histogram.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  -> fig2_histogram.pdf / .png")


def fig3_panels():
    fast_rows = load_csv('helio_panel_fast.csv')
    slow_rows = load_csv('helio_panel_slow.csv')
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5),
                             gridspec_kw={'hspace': 0.48})

    # top panel
    ax = axes[0]
    fast_planets = ['Mercury', 'Venus', 'Earth', 'Mars']
    years_f = np.array([r['year'] for r in fast_rows])

    all_vals = []
    for p in fast_planets:
        all_vals.extend([r[p] for r in fast_rows])
    ypad = 8
    ymin_f = min(all_vals) - ypad
    ymax_f = max(all_vals) + ypad

    for p in fast_planets:
        vals = np.array([r[p] for r in fast_rows])
        mean = vals.mean()
        lw = 1.8 if p == 'Earth' else 0.9
        ls = '--' if p == 'Earth' else '-'
        ax.plot(years_f, vals, color=COLORS[p],
                label=f'{p} (mean {mean:+.1f}' + r'$\degree$)',
                linewidth=lw, linestyle=ls, alpha=0.9)
        ax.axhline(mean, color=COLORS[p], lw=0.5, linestyle=':',
                   alpha=0.5)

    ax.axhline(0, color='black', lw=0.5, linestyle=':')
    ax.set_xlim(0, 5.0)
    ax.set_ylim(ymin_f, ymax_f)
    ax.set_xlabel('Years in series')
    ax.set_ylabel(r'$\delta_k$ (degrees)')
    ax.set_title('Fast planets -- 5-year window (daily resolution)',
                 fontsize=10, pad=4)
    ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.22),
              framealpha=0.95, fontsize=8, borderpad=0.5,
              bbox_transform=ax.transAxes)

    # bottom panel
    ax = axes[1]
    slow_planets = ['Jupiter', 'Saturn']
    years_s = np.array([r['year'] for r in slow_rows])

    all_slow = []
    for p in slow_planets:
        all_slow.extend([r[p] for r in slow_rows])
    ymin_s = min(all_slow) - ypad
    ymax_s = max(all_slow) + 18  # extra headroom for legend above Saturn

    for p in slow_planets:
        vals = np.array([r[p] for r in slow_rows])
        mean = vals.mean()
        ax.plot(years_s, vals, color=COLORS[p], label=p,
                linewidth=1.4, alpha=0.9)
        ax.axhline(mean, color=COLORS[p], lw=0.8, linestyle=':',
                   alpha=0.55,
                   label=f'mean ({mean:+.1f}' + r'$\degree$)')

    ax.axhline(0, color='black', lw=0.5, linestyle=':')

    # annotate Jupiter ~12yr oscillation
    ax.annotate('', xy=(32, ymin_s + 4), xytext=(20, ymin_s + 4),
                arrowprops=dict(arrowstyle='<->', color='#8E44AD', lw=1.0))
    ax.text(26, ymin_s + 2,
            r'$\approx\!12\,\mathrm{yr}$',
            ha='center', fontsize=7.5, color='#8E44AD')

    ax.set_xlim(0, 100)
    ax.set_ylim(ymin_s, ymax_s)
    ax.set_xlabel('Years in series')
    ax.set_ylabel(r'$\delta_k$ (degrees)')
    ax.set_title('Slow planets -- full 100-year series (weekly resolution)',
                 fontsize=10, pad=4)
    ax.legend(ncol=2, loc='upper right', framealpha=0.9)

    fig.text(0.005, 0.5,
             r'$\delta_k(t) = \lambda_k(t) - \lambda_k(t - T^*)$  '
             '(heliocentric ecliptic longitude)',
             va='center', ha='left', rotation='vertical',
             fontsize=7.5, color='#555555')

    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig3_panels.pdf', bbox_inches='tight')
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig3_panels.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  -> fig3_panels.pdf / .png")


def fig4_polar():
    rows    = load_polar()
    planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn']
    epochs  = sorted(set(r['epoch'] for r in rows))

    fig = plt.figure(figsize=(13, 8.5))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.28)

    for ei, epoch in enumerate(epochs):
        ax = fig.add_subplot(gs[ei // 3, ei % 3], projection='polar')
        epoch_rows = {r['planet']: r for r in rows if r['epoch'] == epoch}

        for p in planets:
            pr  = epoch_rows.get(p)
            if pr is None:
                continue
            col = COLORS[p]
            r_p = POLAR_RADIUS[p]

            lon_ref  = np.radians(pr['lon_ref'])
            lon_cand = np.radians(pr['lon_cand'])

            # faint orbit circle
            orb = np.linspace(0, 2 * np.pi, 120)
            ax.plot(orb, [r_p] * 120, color=col, lw=0.25, alpha=0.18,
                    zorder=1)

            # arc between positions — always take the short way round
            diff = lon_ref - lon_cand
            # normalise diff to (-pi, +pi]
            diff = (diff + np.pi) % (2 * np.pi) - np.pi
            n_arc = max(3, int(abs(diff) / np.radians(2)))
            arc = np.linspace(lon_cand, lon_cand + diff, n_arc)
            ax.plot(arc, [r_p] * n_arc, color=col, lw=0.9, alpha=0.5)

            # candidate: open circle
            ax.scatter(lon_cand, r_p, s=40,
                       facecolors='none', edgecolors=col,
                       linewidths=1.4, zorder=4)
            # reference: filled circle
            ax.scatter(lon_ref, r_p, s=55, c=col, zorder=5)

        ax.set_ylim(0, 1.08)
        ax.set_yticks([])
        ax.set_xticks(np.radians([0, 90, 180, 270]))
        ax.set_xticklabels(['0', '90', '180', '270'], fontsize=7)
        ax.tick_params(axis='x', pad=1)
        ax.grid(True, alpha=0.15, linewidth=0.4)

        # planet initials on first panel
        pass  # no orbit labels

        era = 'BCE' if epoch < 0 else 'CE'
        ax.set_title(f'{abs(epoch)} {era}', fontsize=9, pad=7,
                     fontweight='bold')

    # legend cell
    ax_leg = fig.add_subplot(gs[1, 2])
    ax_leg.axis('off')
    from matplotlib.lines import Line2D
    # Main legend: color → planet name only
    handles = []
    for p in planets:
        handles.append(Line2D([0],[0], marker='o', color='w',
                               markerfacecolor=COLORS[p],
                               markeredgecolor=COLORS[p],
                               markersize=9, label=p))
    ax_leg.legend(handles=handles, loc='center', framealpha=0.95,
                  fontsize=9, handletextpad=0.5, borderpad=1.0,
                  title='Planets', title_fontsize=8.5)

    # Caption below figure explaining filled vs open
    fig.text(0.5, -0.02,
             r'Filled circle = heliocentric position at epoch $t$  |  '
             r'Open circle = position at $t - T^*$ ($T^* \approx 1151$ yr)',
             ha='center', fontsize=8.5, color='#444444',
             style='italic')

    fig.suptitle(
        r'Heliocentric positions at epoch $t$ (filled) and $t - T^*$ (open)'
        '\n'
        r'$T^* = 420{,}403\,\mathrm{days} \approx 1151\,\mathrm{yr}$'
        '  --  5 independent epochs spanning 750 years',
        fontsize=9.5, y=1.01
    )

    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig5_polar.pdf', bbox_inches='tight')
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig5_polar.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  -> fig5_polar.pdf / .png")



def fig4_convergence():
    rows = load_csv('helio_convergence.csv')
    lengths = [r['series_years'] for r in rows]
    tstar_scores  = [r['tstar_score']  for r in rows]
    second_scores = [r['second_score'] for r in rows]
    gaps = [s - t for s, t in zip(second_scores, tstar_scores)]

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), gridspec_kw={'wspace': 0.38})

    # Left: T* score and 2nd best vs series length
    ax = axes[0]
    ax.plot(lengths, tstar_scores,  'o-', color='#C0392B', lw=1.4,
            markersize=5, label=r'$T^*$ score')
    ax.plot(lengths, second_scores, 's--', color='#888888', lw=1.2,
            markersize=4, label='2nd best score')
    ax.set_xlabel('Series length (years)')
    ax.set_ylabel(r'Score $S$ (degrees)')
    ax.set_xscale('log')
    ax.set_xticks(lengths)
    ax.set_xticklabels([str(int(l)) for l in lengths], fontsize=8)
    ax.legend(framealpha=0.9)
    ax.set_title('Score convergence vs series length', fontsize=10, pad=4)

    # Right: gap between T* and 2nd best
    ax = axes[1]
    ax.plot(lengths, gaps, 'D-', color='#4A90D9', lw=1.4, markersize=5)
    ax.axhline(np.mean(gaps), color='#4A90D9', lw=0.8, linestyle='--', alpha=0.6,
               label=f'mean gap = {np.mean(gaps):.2f}°')
    ax.set_xlabel('Series length (years)')
    ax.set_ylabel(r'Score gap $S_{2nd} - S_{T^*}$ (degrees)')
    ax.set_xscale('log')
    ax.set_xticks(lengths)
    ax.set_xticklabels([str(int(l)) for l in lengths], fontsize=8)
    ax.legend(framealpha=0.9)
    ax.set_title(r'Gap between $T^*$ and 2nd best candidate', fontsize=10, pad=4)

    fig.suptitle(
        r'$T^* = 420{,}403$ days is the global minimum for all series lengths — '
        r'from 1 to 100 years',
        fontsize=9, y=1.02
    )

    fig.tight_layout()
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig4_convergence.pdf', bbox_inches='tight')
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.0)
    fig.savefig('fig4_convergence.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  -> fig4_convergence.pdf / .png")


def main():
    print("Generating figures for the 1151-year quasi-commensurability paper")
    print("=" * 60)

    required = ['helio_results.csv', 'helio_panel_fast.csv',
                'helio_panel_slow.csv', 'helio_polar.csv',
                'helio_convergence.csv']
    missing = [f for f in required if not os.path.exists(f)]
    if missing:
        print(f"ERROR: Missing input files: {', '.join(missing)}")
        print("Run helio1151.py first to generate these files.")
        return

    print("\n[1/4] Figure 1: Score scatter plot...")
    fig1_scatter()
    print("[2/5] Figure 2: Score histogram...")
    fig2_histogram()
    print("[3/5] Figure 3: Per-planet offset panels...")
    fig3_panels()
    print("[4/5] Figure 4: Polar snapshots...")
    fig4_polar()

    print("[5/5] Figure 5: Series length convergence...")
    if os.path.exists('helio_convergence.csv'):
        fig4_convergence()
    else:
        print("  helio_convergence.csv not found — skipping.")

    print("\n" + "=" * 60)
    print("All figures saved as PDF and PNG.")
    print("Upload the PDF files to Overleaf.")
    print("=" * 60)


if __name__ == '__main__':
    main()
