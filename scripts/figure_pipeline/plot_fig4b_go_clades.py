#!/usr/bin/env python3
"""plot_fig4b_go_clades.py - GO enrichment of the corymbose-specific,
solitary-specific and shared (corymbose & solitary) family sets, using the
all-17-track annotation.

Bubble area scales with sqrt(number of families carrying the term); the label
gives families/background; a size legend is drawn below each panel.

Inputs (analysis/racemose_inflorescence_scan_20260922/results_v2/):
  go_enrichment_corymbose_all17.tsv
  go_enrichment_solitary_all17.tsv
  go_enrichment_cor_sol_shared_all17.tsv

Output: fig4b_go_cor_sol_shared.{png,pdf,svg}
"""
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

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
RES = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'

SIZE_BASE, SIZE_SCALE = 12.0, 30.0


def bubble_size(k):
    return SIZE_BASE + SIZE_SCALE * math.sqrt(k)


def size_legend_handles(color, ks=(5, 50, 200)):
    return [Line2D([], [], marker='o', ls='', color=color, alpha=0.85,
                   ms=2 * math.sqrt(bubble_size(k) / math.pi), label=str(k))
            for k in ks]


def go_bubble_panel(ax, rows, color, label_len=48, xlabel='\u2212log10 FDR (BH)'):
    y = np.arange(len(rows))
    x = np.array([-math.log10(max(float(r['q']), 1e-12)) for r in rows])
    k = np.array([int(r['expanded']) for r in rows])
    s = np.array([bubble_size(v) for v in k])
    ax.scatter(x, y, s=s, color=color, alpha=0.85, lw=0.5, edgecolor='white', zorder=3)
    xr = max(x.max() - x.min(), 0.5)
    lo = max(0.0, x.min() - 0.10 * xr)
    hi = x.max() + 0.30 * xr
    ax.set_xlim(lo, hi)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    for yy, r in zip(y, rows):
        ax.text(hi - 0.015 * (hi - lo), yy, f"{r['expanded']}/{r['bg']}",
                fontsize=5.4, va='center', ha='right', color='#555555')
    ax.set_yticks(y, [r['go_description'][:label_len] for r in rows], fontsize=5.8)
    ax.set_xlabel(xlabel, fontsize=7.0)
    ax.tick_params(labelsize=6.2)
    ax.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    return ax


def load(path, top, qmax=0.05):
    rows = [r for r in csv.DictReader(path.open(), delimiter='\t')]
    rows = [r for r in rows if float(r['q']) < qmax][:top]
    return rows[::-1]


def main():
    panels = [
        ('a  Corymbose-specific (1,255 families)', 'go_enrichment_corymbose_all17.tsv', 10, '#E15759'),
        ('b  Solitary-specific (1,054 families)', 'go_enrichment_solitary_all17.tsv', 10, '#4E79A7'),
        ('c  Shared corymbose & solitary (80 families)', 'go_enrichment_cor_sol_shared_all17.tsv', 8, '#59A14F'),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(5.2, 6.4),
                             gridspec_kw=dict(hspace=0.62))
    for i, (ax, (title, fname, top, color)) in enumerate(zip(axes, panels)):
        rows = load(RES / fname, top)
        go_bubble_panel(ax, rows, color)
        ax.set_title(title, fontsize=8.0, loc='left')
    handles = size_legend_handles('#8C8C8C')
    fig.legend(handles=handles, title='families with term', fontsize=5.4,
               title_fontsize=5.6, loc='lower right', bbox_to_anchor=(1.0, -0.02),
               ncol=3, columnspacing=2.2, handletextpad=0.5, handleheight=2.4)
    fig.text(0.01, 0.005,
             'bubble area \u221d \u221a(families with term); label = families/background; '
             'dotted line: FDR = 0.05',
             fontsize=5.2, color='#555555')
    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        fig.savefig(FIG / f'fig4b_go_cor_sol_shared.{ext}', bbox_inches='tight', **kw)
    plt.close(fig)
    print('fig4b_go_cor_sol_shared saved')


if __name__ == '__main__':
    main()
