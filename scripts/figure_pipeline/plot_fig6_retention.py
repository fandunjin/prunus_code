#!/usr/bin/env python3
"""plot_fig6_retention.py - manuscript Figure 6 for the post-WGD retention result
group, plus the two former fig7 panels drawn separately (no distortion).

fig6_retention panels
  a  CAFE5 family-size change along the dated tree        plot_figures_v2.fig1
  b  genome-wide scan volcano                             plot_fig2c_scan_all17.fig2c
  c  consistent racemose-retained families (heatmap)      plot_all17_figures.fig3_all17
  d  GO enrichment of the retained set                    plot_all17_figures.fig4_all17
  e  set overlap of the three clade sets (Venn, round)    this file
  f  branch-model dN/dS significance rate                 this file

Separate panels (former fig7, now split and aspect-corrected):
  fig7_shared_specific_all17a.{png,pdf,svg}   Venn
  fig7_shared_specific_all17b.{png,pdf,svg}   GO profiles of the clade sets
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import plot_figures_v2 as pv                     # noqa: E402
import plot_all17_figures as paf                 # noqa: E402
from plot_fig2c_scan_all17 import fig2c, fig2c_broken   # noqa: E402
from plot_fig4b_go_clades import (go_bubble_panel, size_legend_handles)  # noqa: E402
from plot_fig4_go_grouped import go_grouped_panel                        # noqa: E402

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

_CANDIDATES = [Path(r'${PROJECT_ROOT}'),
               Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'
RES = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'
DNDS = B / 'analysis/dnds_branch_20260922/results'


def save(fig, name):
    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        fig.savefig(FIG / f'{name}.{ext}', bbox_inches='tight', **kw)
    plt.close(fig)
    print('saved', name)


# ------------------------------------------------------------------ panels
def venn_panel(ax, title=None):
    """Overlap of the three clade sets; equal aspect keeps the circles round."""
    counts = {}
    with (RES / 'shared_specific_counts_all17.tsv').open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['category']] = int(r['n_families'])
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.45, 1.25)
    ax.set_aspect('equal')
    ax.axis('off')
    for (cx, cy), col in [((-0.42, 0.28), '#C0504D'), ((0.42, 0.28), '#4F81BD'),
                          ((0.0, -0.42), '#59A14F')]:
        ax.add_patch(Circle((cx, cy), 0.78, fc=col, alpha=0.20, ec=col, lw=1.0))
    put = lambda x, y, n, fs=7.2: ax.text(x, y, str(n), ha='center', va='center',
                                          fontsize=fs)
    put(-0.92, 0.32, counts['R_only'])
    put(0.92, 0.32, counts['C_only'])
    put(0.0, -0.92, counts['S_only'])
    put(0.0, 0.58, counts['R_C'], 6.4)
    put(-0.60, -0.50, counts['R_S'], 6.4)
    put(0.60, -0.50, counts['C_S'], 6.4)
    put(0.0, -0.16, counts['R_C_S'], 6.4)
    ax.text(-1.08, 0.94, f'racemose\n(n={counts["R_only"] + counts["R_C"] + counts["R_S"] + counts["R_C_S"]:,})', color='#C0504D', fontsize=6.4,
            ha='center')
    ax.text(1.08, 0.94, f'corymbose\n(n={counts["C_only"] + counts["R_C"] + counts["C_S"] + counts["R_C_S"]:,})', color='#4F81BD', fontsize=6.4,
            ha='center')
    ax.text(0.0, -1.30, f'solitary (n={counts["S_only"] + counts["R_S"] + counts["C_S"] + counts["R_C_S"]:,})', color='#59A14F', fontsize=6.4,
            ha='center')
    if title:
        ax.set_title(title, fontsize=8.4, loc='left')


def clade_go_panel(ax, title=None, label_len=38,
                   picks=(('R_only', 4), ('C_only', 4), ('S_only', 4), ('C_S', 3))):
    """GO profiles of the clade sets (former fig7 panel b)."""
    rows = list(csv.DictReader((RES / 'shared_specific_go_all17.tsv').open(),
                               delimiter='\t'))
    cat_col = {'R_only': '#C0504D', 'C_only': '#4F81BD', 'S_only': '#59A14F',
               'C_S': '#9A4D8E'}
    sel = []
    for g, n in picks:
        rs = sorted([r for r in rows if r['group'] == g],
                    key=lambda r: float(r['q']))[:n]
        sel += rs[::-1]
    ys = np.arange(len(sel))
    for y, r in zip(ys, sel):
        col = cat_col[r['group']]
        x = -math.log10(max(float(r['q']), 1e-12))
        ax.plot([0, x], [y, y], color=col, lw=0.8, alpha=0.55)
        ax.scatter([x], [y], s=10 + 4 * int(r['k']), color=col, zorder=3)
    ax.set_yticks(ys, [r['go_description'][:label_len] for r in sel], fontsize=5.2)
    for tick, r in zip(ax.get_yticklabels(), sel):
        tick.set_color(cat_col[r['group']])
    ax.set_xlabel('\u2212log10 q (BH)', fontsize=7.0)
    ax.tick_params(axis='x', labelsize=6.4)
    ax.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    ax.set_xlim(0, max(-math.log10(max(float(r['q']), 1e-12)) for r in sel) * 1.12)
    handles = [Line2D([], [], marker='o', ls='', color=c, label=l, ms=4)
               for l, c in [('racemose-only', '#C0504D'), ('corymbose-only', '#4F81BD'),
                            ('solitary-only', '#59A14F'), ('shared cor & sol', '#9A4D8E')]]
    ax.legend(handles=handles, fontsize=5.4, loc='lower right')
    if title:
        ax.set_title(title, fontsize=8.4, loc='left')


def dnds_panel(ax, title=None):
    """Branch-model dN/dS: fraction of sets with p < 0.05, racemose vs background."""
    groups = defaultdict(list)
    for r in csv.DictReader((DNDS / 'branch_dnds_integrated.tsv').open(),
                            delimiter='\t'):
        groups[r['group']].append(r)
    order = [('background_A', 'background, A copies', '#A8A8A8'),
             ('background_B', 'background, B copies', '#D8D8D8'),
             ('control', 'family controls', '#4F81BD'),
             ('candidate', 'candidates', '#C0504D')]
    ys = np.arange(len(order))
    for y, (g, lab, col) in zip(ys, order):
        rs = groups[g]
        frac = sum(1 for r in rs if float(r['p']) < 0.05) / len(rs)
        wr = np.mean([float(r['omega_rac']) for r in rs])
        wb = np.mean([float(r['omega_bg']) for r in rs])
        ax.barh(y, frac, height=0.62, color=col, zorder=3)
        ax.text(frac + 0.012, y, f'{frac * 100:.0f}%   (\u03c9 {wr:.2f} vs {wb:.2f})',
                va='center', fontsize=5.6, color='#3C3C3C')
    ax.axvline(0.05, color='#C0504D', ls=':', lw=0.9, zorder=4)
    ax.text(0.058, 3.46, '5% null', fontsize=5.2, color='#C0504D', va='center')
    ax.set_yticks(ys, [lab for _, lab, _ in order], fontsize=6.0)
    ax.set_xlim(0, 0.62)
    ax.set_ylim(-0.6, 3.62)
    ax.set_xlabel('fraction of sets with p < 0.05', fontsize=7.0)
    ax.tick_params(labelsize=6.2)
    if title:
        ax.set_title(title, fontsize=8.4, loc='left')


# ------------------------------------------------------------------ figures
def fig7_split():
    fig = plt.figure(figsize=(3.5, 3.3))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    venn_panel(ax, 'a  Overlap of the enrichment sets')
    save(fig, 'fig7_shared_specific_all17a')

    fig = plt.figure(figsize=(5.2, 3.6))
    ax = fig.add_axes([0.33, 0.13, 0.64, 0.80])
    clade_go_panel(ax, 'b  GO profiles of the clade sets (all-17 annotation)',
                   label_len=44)
    fig.text(0.02, 0.02, 'dot size \u221d families with term', fontsize=5.2,
             color='#555555')
    save(fig, 'fig7_shared_specific_all17b')


def fig6_retention():
    fig = plt.figure(figsize=(178 / 25.4, 182 / 25.4))
    gs = fig.add_gridspec(4, 3, left=0.165, right=0.988, top=0.965, bottom=0.075,
                          hspace=0.55, wspace=0.42,
                          width_ratios=[1.0, 1.0, 0.74],
                          height_ratios=[1.0, 1.0, 0.80, 0.52])

    axa = fig.add_subplot(gs[0, 0:2])
    pv.fig1(axa, compact=True)
    axa.set_title('a  CAFE5 family-size change along the dated tree', fontsize=8.2,
                  loc='left')

    axb = fig.add_subplot(gs[0, 2])
    venn_panel(axb, 'b  Set overlap')

    sub = gs[1, 0:2].subgridspec(2, 1, height_ratios=[0.45, 2.45], hspace=0.0)
    axc_hi = fig.add_subplot(sub[0, 0])
    axc = fig.add_subplot(sub[1, 0], sharex=axc_hi)
    h, l = fig2c_broken(axc, axc_hi, s=8.0, annotate=True, legend=True,
                        yticks_lo=(0, 20, 40), yticks_hi=(60, 120))
    axc.set_ylabel('\u2212log10 p', fontsize=6.6)
    axc.legend(h, l, fontsize=4.0, loc='upper left', ncol=2, columnspacing=0.6,
               handletextpad=0.2, labelspacing=0.22)
    axc_hi.set_title('c  Genome-wide scan (all-17 annotation)', fontsize=8.2,
                     loc='left')

    axd = fig.add_subplot(gs[1, 2])
    go_grouped_panel(axd, size_base=4.0, size_scale=7.0, label_fs=5.2,
                     show_legend=True, n_per_group=2, show_counts=False)
    axd.set_title('d  GO enrichment', fontsize=8.2, loc='left')

    axe = fig.add_subplot(gs[2, 0:3])
    paf.fig3_all17(axe)
    axe.set_title('e  Consistent racemose-retained families (copies, enrichment, '
                  'significance)', fontsize=8.2, loc='left', pad=16)

    axf = fig.add_subplot(gs[3, 0:3])
    dnds_panel(axf, 'f  Branch-model dN/dS on the racemose branches')

    fig.text(0.165, 0.006,
             'a: +up (\u2212down) families per branch; parentheses = branch-significant '
             '(p < 0.05).\n'
             'b: set sizes and overlaps of the three clade sets; c: log2 enrichment '
             'vs corymbose + solitary.\n'
             'd: bubble area \u221d \u221a(families with term), label = families/background.\n'
             'e: right panels = enrichment and significance; stars = 2/3 (**) or 3/3 '
             '(***) racemose species; f: LRT vs single-\u03c9 null.',
             fontsize=4.2, color='#555555', va='bottom', linespacing=1.5)
    save(fig, 'fig6_retention')


if __name__ == '__main__':
    fig7_split()
    fig6_retention()
