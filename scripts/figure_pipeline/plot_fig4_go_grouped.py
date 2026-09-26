#!/usr/bin/env python3
"""plot_fig4_go_grouped.py - GO enrichment of the racemose-retained families,
grouped by GO namespace (biological process, cellular component, molecular
function).

The three namespaces are separated by colour (bubbles + term labels) and by
thin separators; the bubble area scales with sqrt(families with term) and the
label gives families/background.

Input : analysis/racemose_inflorescence_scan_20260922/results_v2/
        go_enrichment_racemose_all17.tsv
        0903/annotation/functional/pzi/swissprot_go.tsv   (go_id -> aspect)
Output: fig4_go_enrichment_all17.{png,pdf,svg}

go_grouped_panel(ax, ...) is reused by plot_fig6_retention.py (Fig. 6d).
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
    'axes.linewidth': 0.6,
    'legend.frameon': False,
})

_CANDIDATES = [Path(r'${PROJECT_ROOT}'),
               Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'
RES = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'
GO_SRC = B / '0903/annotation/functional/pzi/swissprot_go.tsv'

C_BP = '#C0504D'      # biological process
C_CC = '#4F81BD'      # cellular component
C_MF = '#59A14F'      # molecular function

ASPECTS = [('P', 'biological process', C_BP),
           ('C', 'cellular component', C_CC),
           ('F', 'molecular function', C_MF)]

# concise display labels (full names stay in the source table)
SHORT = {
    'DNA-binding transcription repressor activity, RNA polymerase II-specific':
        'transcription repressor',
    'endoplasmic reticulum to Golgi vesicle-mediated transport':
        'ER\u2192Golgi transport',
    'regulation of alternative mRNA splicing, via spliceosome':
        'mRNA splicing (alt.)',
    'homologous chromosome pairing at meiosis': 'chromosome pairing',
    'meiotic DNA double-strand break formation': 'meiotic DSB formation',
    'vegetative to reproductive phase transition of meristem':
        'vegetative\u2192reproductive',
    'mRNA splicing, via spliceosome': 'mRNA splicing',
    'ATP-dependent protein folding chaperone': 'protein folding chaperone',
    'primary miRNA processing': 'pri-miRNA processing',
    'nuclear dicing body': 'nuclear dicing body',
    'spliceosomal complex': 'spliceosomal complex',
    'Golgi trans cisterna': 'Golgi trans cisterna',
    'trans-Golgi network': 'trans-Golgi network',
}


def load_terms(n_per_group=None):
    aspect = {}
    with GO_SRC.open(encoding='utf-8', errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            aspect.setdefault(r['go_id'], r['go_aspect'])
    rows = list(csv.DictReader((RES / 'go_enrichment_racemose_all17.tsv').open(),
                               delimiter='\t'))
    out = []
    for code, name, color in ASPECTS:
        sel = sorted([r for r in rows if aspect.get(r['go_id']) == code],
                     key=lambda r: float(r['q']))
        if n_per_group is not None:
            sel = sel[:n_per_group]
        out.append((name, color, sel))
    return out


def go_grouped_panel(ax, size_base=6.0, size_scale=15.0, label_fs=5.6,
                     show_legend=False, gap=0.7, n_per_group=None,
                     show_counts=True, count_fs=None):
    """Grouped GO bubble panel (three GO namespaces)."""
    groups = load_terms(n_per_group)
    ys, marks = [], []
    y = 0.0
    for name, color, terms in groups:
        first = y
        for r in terms:
            ys.append((y, r, color))
            y += 1.0
        marks.append((name, color, first - 0.35, y - 0.65))
        y += gap
    yy = np.array([p[0] for p in ys])
    x = np.array([-math.log10(max(float(r['q']), 1e-12)) for _, r, _ in ys])
    k = np.array([int(r['expanded']) for _, r, _ in ys])
    s = np.array([size_base + size_scale * math.sqrt(v) for v in k])

    lo = min(-math.log10(0.05), x.min()) - 0.12
    hi = x.max() + (0.10 if show_counts else 0.06)
    for (y0, r, color), xi, si in zip(ys, x, s):
        ax.plot([lo, xi], [y0, y0], color=color, lw=0.5, alpha=0.35, zorder=2)
        ax.scatter([xi], [y0], s=si, color=color, alpha=0.9, lw=0.5,
                   edgecolor='white', zorder=3)
        if show_counts:
            ax.text(hi + 0.04, y0, f"{r['expanded']}/{r['bg']}",
                    fontsize=count_fs or (label_fs - 0.6), va='center',
                    ha='left', color='#666666', clip_on=False)
    ax.set_yticks(yy, [SHORT.get(r['go_description'], r['go_description'])
                       for _, r, _ in ys], fontsize=label_fs)
    for tick, (_, _, color) in zip(ax.get_yticklabels(), ys):
        tick.set_color(color)
    ax.set_ylim(y - gap + 0.4, -0.8)
    ax.set_xlim(lo, hi)
    ax.axvline(-math.log10(0.05), color='#999999', ls=':', lw=0.6, zorder=1)
    ax.set_xlabel('\u2212log10 FDR (BH)', fontsize=7.0, labelpad=2)
    ax.tick_params(axis='x', labelsize=6.2, length=2.0, width=0.6)
    ax.tick_params(axis='y', length=0)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(0.6)

    # thin separators between the three namespaces
    for (name, color, y0, y1), (n2, c2, y2, y3) in zip(marks[:-1], marks[1:]):
        ax.axhline((y1 + y2) / 2.0, color='#DDDDD8', lw=0.6, zorder=1)

    if show_legend:
        handles = [Line2D([], [], marker='o', ls='', color=c, ms=4, label=n)
                   for n, c, _, _ in marks]
        ax.legend(handles=handles, fontsize=5.0, loc='lower right',
                  bbox_to_anchor=(0.995, 0.02), handletextpad=0.4,
                  labelspacing=0.35, borderpad=0.2)
    return marks


SB, SS = 5.5, 13.0          # bubble size scale (area = SB + SS*sqrt(k))


def main():
    fig = plt.figure(figsize=(3.10, 4.80))
    ax = fig.add_axes([0.445, 0.265, 0.495, 0.675])
    marks = go_grouped_panel(ax, size_base=SB, size_scale=SS, label_fs=5.0,
                             n_per_group=6, count_fs=4.6)

    group_handles = [Line2D([], [], marker='o', ls='', color=c, ms=4.0,
                            label=lab) for (n, c, _, _), lab in
                     zip(marks, ('BP', 'CC', 'MF'))]
    fig.legend(handles=group_handles, fontsize=4.8, loc='lower left',
               bbox_to_anchor=(0.445, 0.150), ncol=3, columnspacing=1.0,
               handletextpad=0.3)
    size_handles = [Line2D([], [], marker='o', ls='', color='#8C8C8C', alpha=0.85,
                           ms=2 * math.sqrt((SB + SS * math.sqrt(k)) / math.pi),
                           label=str(k)) for k in (5, 50, 200)]
    ax.legend(handles=size_handles, fontsize=4.4, loc='lower right',
              bbox_to_anchor=(0.995, 0.012), ncol=3, columnspacing=0.8,
              handletextpad=0.25, borderpad=0.15, labelspacing=0.4)
    fig.text(0.445, 0.105,
             'BP, biological process; CC, cellular component;',
             fontsize=4.0, color='#555555')
    fig.text(0.445, 0.078, 'MF, molecular function.', fontsize=4.0,
             color='#555555')
    fig.text(0.445, 0.018,
             'bubble area \u221d \u221a(families with term);\n'
             'label = families/background; dotted line: FDR = 0.05.',
             fontsize=4.0, color='#555555', linespacing=1.45)
    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        fig.savefig(FIG / f'fig4_go_enrichment_all17.{ext}', bbox_inches='tight', **kw)
    plt.close(fig)
    print('saved fig4_go_enrichment_all17')


if __name__ == '__main__':
    main()
