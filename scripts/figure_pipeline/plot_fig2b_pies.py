#!/usr/bin/env python3
"""plot_fig2b_pies.py - volcano with copy-source pie markers (fig2b).

Same statistics as fig2_scan_volcano (plot_figures_v2.py::load_scan):
  x = log2 enrichment (racemose vs corymbose+solitary tracks)
  y = -log10 min(p_greater, p_less), binomial test, n = 186,830 racemose copies
Background = all tested families (small gray dots).
Foreground = the top `--top-n` families by significance drawn as pie markers,
slices = copy source (PZI / PSP / PPA / corymbose / solitary / outgroup),
pie area proportional to the family's total copy number.
Broken y-axis (0-36 | 36-120) so the dense significance tail stays readable.

Outputs (0903/figures/inflorescence_evolution_review/):
  fig2b_scan_volcano_pies.{png,pdf,svg}
  data/fig2b_pie_selection.tsv   (source data for the pie selection)
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.patches import Circle, Wedge

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'font.size': 7,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.linewidth': 0.8,
    'legend.frameon': False,
})

_CANDIDATES = [Path(r'${PROJECT_ROOT}'), Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'

TOT = {'RAC': 186830, 'COR': 131314, 'SOL': 150034}

GROUPS = [
    ('PZI', 'pzi', '#B64342'),
    ('PSP', 'psp', '#E28E2C'),
    ('PPA', 'ppa', '#9A4D8E'),
    ('Corymbose', 'cor', '#5B8FD6'),
    ('Solitary', 'sol', '#7BAA5B'),
    ('Outgroup', 'out', '#8F8F8F'),
]

HIGHLIGHT = {'OG0000961': 'BLH1', 'OG0005156': 'ASY1', 'OG0005416': 'PDR2',
             'OG0000864': 'CPK1', 'OG0000130': 'LNK1', 'OG0000222': 'RAPTOR2',
             'OG0001188': 'Susy', 'OG0000664': 'enolase', 'OG0000348': 'MCM7',
             'OG0000256': 'pol (TE)'}

OFFSETS = {
    'OG0000256': (-6, 2), 'OG0000130': (-6, 2), 'OG0000228': (-7, 2),
    'OG0001439': (5, 2), 'OG0000222': (-6, 3), 'OG0000664': (1, 6),
    'OG0000864': (8, 5), 'OG0000348': (-9, 5), 'OG0003061': (5, -9),
    'OG0000515': (-9, -8), 'OG0000091': (6, -8), 'OG0000129': (-7, 4),
    'OG0000961': (-9, -9), 'OG0005416': (4, 3), 'OG0001188': (-8, -7),
    'OG0005156': (6, -9),
}


def binom_pmf_log(n, k, p):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def binom_less(k, n, p):
    if k >= n or p >= 1:
        return 1.0
    if p <= 0:
        return 0.0
    return sum(math.exp(binom_pmf_log(n, i, p)) for i in range(0, k + 1))


def load_points():
    rows = list(csv.DictReader((FIG / 'data/clade_family_scan_v3.tsv').open(), delimiter='\t'))
    pts = []
    for r in rows:
        if float(r['rac']) + float(r['cor']) + float(r['sol']) + float(r['out']) == 0:
            continue
        try:
            enr = float(r['enr_cs'])
        except Exception:
            continue
        if not math.isfinite(enr) or enr <= 0:
            continue
        base = float(r['cor']) + float(r['sol'])
        frac = base / (TOT['COR'] + TOT['SOL'])
        p_up = float(r['p_cs'])
        if enr >= 1.0 or base <= 0:
            p_dn = 1.0
        else:
            p_dn = binom_less(int(r['rac']), TOT['RAC'], frac)
        p = min(p_up, p_dn)
        total = float(r['rac']) + float(r['cor']) + float(r['sol']) + float(r['out'])
        pts.append(dict(r, x=math.log2(enr), y=-math.log10(max(p, 1e-300)),
                        p=p, total=total))
    return pts


def add_pie(ax, xy, fracs, r_pt, zorder, xycoords='data', te=False, lw_edge=0.35,
            flat_color=None):
    da = DrawingArea(2 * r_pt, 2 * r_pt, 0, 0)
    start = 90.0
    for (name, key, color), f in zip(GROUPS, fracs):
        if f <= 0:
            continue
        dtheta = 360.0 * f
        da.add_artist(Wedge((r_pt, r_pt), r_pt, start - dtheta, start,
                            facecolor=(flat_color or color),
                            edgecolor='white', linewidth=lw_edge))
        start -= dtheta
    if te:
        da.add_artist(Circle((r_pt, r_pt), r_pt, facecolor='none', edgecolor='#4D4D4D',
                             linewidth=0.7, linestyle=(0, (1.4, 0.9))))
    ab = AnnotationBbox(da, xy, xycoords=xycoords, frameon=False, pad=0,
                        box_alignment=(0.5, 0.5))
    ab.set_zorder(zorder)
    ax.add_artist(ab)


def fig2b(top_n=60, pie_scale=1.0, split=36.0, figsize=(7.4, 5.4)):
    pts = load_points()
    sel = sorted(pts, key=lambda d: d['y'], reverse=True)[:top_n]

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 4.0], hspace=0.06,
                          left=0.075, right=0.985, top=0.925, bottom=0.10)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
    ax_bot.set_ylim(0, split)
    ax_top.set_ylim(split, 120)
    ax_top.spines['bottom'].set_visible(False)
    ax_top.tick_params(bottom=False, labelbottom=False)
    ax_bot.set_yticks([0, 10, 20, 30])
    ax_top.set_yticks([40, 60, 80, 100])

    d = 0.012
    for ax, edge in ((ax_top, 0.0), (ax_bot, 1.0)):
        ax.plot((-d, +d), (edge - d * 2, edge + d * 2), transform=ax.transAxes,
                color='black', clip_on=False, lw=0.8)

    for ax in (ax_top, ax_bot):
        ax.scatter([p['x'] for p in pts], [p['y'] for p in pts], s=1.1,
                   c='#D9D9D9', alpha=0.55, lw=0, zorder=1)
        ax.axvline(0, color='#888888', lw=0.6, ls='--', zorder=2)
    ax_bot.axhline(-math.log10(0.05), color='#888888', lw=0.6, ls=':', zorder=2)

    t_lo, t_hi = 3.0, 150.0
    r_lo, r_hi = 3.0 * pie_scale, 6.8 * pie_scale

    def radius(total):
        f = (math.sqrt(max(total, t_lo)) - math.sqrt(t_lo)) / (math.sqrt(t_hi) - math.sqrt(t_lo))
        return r_lo + (r_hi - r_lo) * min(max(f, 0.0), 1.0)

    draw = sorted(sel, key=lambda dd: (dd['y'], dd['total']))
    for i, dd in enumerate(draw):
        counts = [float(dd[k]) for _, k, _ in GROUPS]
        s = sum(counts)
        fracs = [c / s if s > 0 else 0.0 for c in counts]
        ax = ax_top if dd['y'] > split else ax_bot
        add_pie(ax, (dd['x'], dd['y']), fracs, radius(dd['total']),
                zorder=5 + i * 0.01, te=(dd['te17'] == 'True'))

    # labels
    labelled = set()
    for dd in sorted(sel, key=lambda x: x['y']):
        og = dd['og']
        text = HIGHLIGHT.get(og)
        if text is None and dd['te17'] != 'True' and (dd['symbol17'] or '').strip() and dd['y'] >= 20:
            text = dd['symbol17'].strip().rstrip(';')
        if not text or og in labelled:
            continue
        labelled.add(og)
        dx, dy = OFFSETS.get(og, (5, 3))
        ha = 'left' if dx > 0 else 'right'
        ax = ax_top if dd['y'] > split else ax_bot
        ax.annotate(text, (dd['x'], dd['y']), xytext=(dx, dy),
                    textcoords='offset points', ha=ha, va='bottom',
                    fontsize=5.6, color='#333333', zorder=40,
                    bbox=dict(fc='white', ec='none', alpha=0.65, pad=0.12),
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='#999999',
                                    shrinkA=0, shrinkB=1))

    # axes cosmetics
    ax_bot.set_xlabel('log2 enrichment (racemose vs diploid Prunus)', fontsize=7.2)
    ax_bot.tick_params(labelsize=6.4)
    ax_top.tick_params(labelsize=6.4)
    fig.supylabel('\u2212log10 p', fontsize=7.2, x=0.012)
    ax_top.set_title(f'b  Genome-wide scan: top {top_n} families by significance '
                     '(pie = copy source)', fontsize=8.4, loc='left')

    # colour legend + size legend (bottom panel, empty upper-left region)
    handles = [Line2D([], [], marker='o', ls='', ms=4.4, color=c, label=n)
               for n, _, c in GROUPS]
    ax_bot.legend(handles=handles, fontsize=5.4, loc='upper left',
                  bbox_to_anchor=(0.005, 0.995), ncol=2, columnspacing=0.8,
                  handletextpad=0.25, labelspacing=0.35)
    ex_x, ex_y0, dy = 0.075, 0.72, 0.075
    ax_bot.text(ex_x, ex_y0 + 0.075, 'family size', transform=ax_bot.transAxes,
                fontsize=5.2, va='bottom', color='#666666')
    for j, total in enumerate((16, 60, 150)):
        yy = ex_y0 - j * dy
        add_pie(ax_bot, (ex_x, yy), [1, 0, 0, 0, 0, 0], radius(total), zorder=30,
                xycoords='axes fraction', flat_color='#9E9E9E')
        ax_bot.text(ex_x + 0.033, yy, f'{total} copies', transform=ax_bot.transAxes,
                    fontsize=5.2, va='center', color='#666666')
    ax_bot.text(0.42, 0.99,
                f'background: all {len(pts):,} tested families\n'
                'pie area scales with total copies\n'
                'dashed ring = TE-derived; colours as legend',
                transform=ax_bot.transAxes, fontsize=5.2, color='#666666', va='top')

    # source data
    out_tsv = FIG / 'data/fig2b_pie_selection.tsv'
    with out_tsv.open('w', newline='', encoding='utf-8') as fh:
        cols = ['rank', 'og', 'symbol17', 'te17', 'pzi', 'psp', 'ppa', 'cor', 'sol', 'out',
                'total', 'x_log2enr', 'y_neglog10p', 'enr_cs', 'p_cs', 'n_species_up']
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for rank, dd in enumerate(sel, 1):
            w.writerow([rank, dd['og'], dd['symbol17'], dd['te17'], dd['pzi'], dd['psp'],
                        dd['ppa'], dd['cor'], dd['sol'], dd['out'], int(dd['total']),
                        f"{dd['x']:.4g}", f"{dd['y']:.4g}", dd['enr_cs'], dd['p_cs'],
                        dd['n_species_up']])

    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        fig.savefig(FIG / f'fig2b_scan_volcano_pies.{ext}', bbox_inches='tight', **kw)
    plt.close(fig)
    print(f'saved fig2b_scan_volcano_pies.{{png,pdf,svg}}; selection y range '
          f'{sel[-1]["y"]:.2f}-{sel[0]["y"]:.2f}; source data -> {out_tsv.name}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--top-n', type=int, default=60)
    ap.add_argument('--pie-scale', type=float, default=1.0)
    args = ap.parse_args()
    fig2b(top_n=args.top_n, pie_scale=args.pie_scale)
