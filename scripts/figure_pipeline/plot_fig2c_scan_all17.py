#!/usr/bin/env python3
"""plot_fig2c_scan_all17.py - genome-wide scan volcano coloured by functional
category using the all-17-track annotation (clade_family_scan_v3.tsv).

Same statistics as fig2_scan_volcano:
  x = log2 enrichment (racemose vs corymbose+solitary tracks)
  y = -log10 min(p_greater, p_less), binomial test, n = 186,830 racemose copies

Output: fig2c_scan_volcano_all17.{png,pdf,svg}
"""
import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

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

CATS = [
    ('TE / transposon', '#E8A33D',
     r'transposon|retrotransposon|retrovirus|transposase|integrase|reverse transcriptase|'
     r'pol polyprotein|gag|gypsy|copia|rnase h|retro'),
    ('floral / meristem / hormone', '#C0504D',
     r'flower|floral|infloresc|meristem|mads|knox|bel1|bell|svp|short vegetative|dam|'
     r'spl9|yabby|ap1|ap2|leafy|wus|clavata|tfl|ft-like|sepallata|soc1|agamous|pistillata|'
     r'auxin|cytokinin|gibberellin|abscisic|ethylene|jasmon|brassinosteroid|strigolactone|'
     r'photoperiod|circadian|dormancy|vernal|flowering'),
    ('defense / stress', '#8C564B',
     r'defense|disease|pathogen|salicylic|jasmonic|hypersensitive|'
     r'heat shock|cold|freez|drought|salt stress|oxidative stress|wound|'
     r'nbs-lrr|nb-arc|tmv resistance|tir domain|receptor-like protein kinase'),
    ('RNA / splicing', '#4F81BD',
     r'splic|rna helicase|deah|dead-box|mrna|sr protein|sc35|ribonucleoprotein|'
     r'polyadenyl|cleavage factor|exon junction|pentatricopeptide'),
    ('transport / vesicle', '#59A14F',
     r'golgi|vesicle|sec13|sec23|coatomer|snare|vacuolar|transporter|permease|'
     r'potassium|sulfate|oligopeptide|er to golgi|emp24'),
    ('meiosis / DNA / chromatin', '#B07AA1',
     r'meiosis|chiasma|mcm|replication|pif1|dna polymerase|recombination|rad51|'
     r'synaptonemal|asy1|dna repair|histone|chromatin|topoisomerase|helicase'),
    ('metabolism / other', '#8C8C8C',
     r'enolase|sucrose synthase|phosphatase|kinase|oxidase|reductase|synthase|transferase|'
     r'protease|lipase|dehydrogenase|peroxidase|methyltransferase|ubiquitin|f-box'),
]

HIGHLIGHT = {'OG0000961': 'BLH1', 'OG0005156': 'ASY1', 'OG0005416': 'PDR2',
             'OG0000864': 'CPK1', 'OG0000130': 'LNK1', 'OG0000222': 'RAPTOR2',
             'OG0001188': 'Susy', 'OG0000664': 'enolase', 'OG0000348': 'MCM7',
             'OG0000256': 'pol (TE)'}

OFFSETS = {
    'OG0000256': (-6, 2), 'OG0000130': (-6, 2), 'OG0000222': (-6, 3),
    'OG0000664': (1, 6), 'OG0000864': (8, 5), 'OG0000348': (-9, 5),
    'OG0000961': (-9, -9), 'OG0005416': (4, 3), 'OG0001188': (-8, -7),
    'OG0005156': (6, -9),
}


def classify(sym, desc):
    text = f'{sym} {desc}'.lower()
    for name, color, pat in CATS:
        if re.search(pat, text):
            return name, color
    return ('annotated (other)', '#BBBBBB') if (sym or desc) else ('unannotated', '#E6E6E6')


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
    out = []
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
        p_dn = binom_less(int(r['rac']), TOT['RAC'], frac) if (enr < 1 and base > 0) else 1.0
        p = min(p_up, p_dn)
        cat, color = classify(r['symbol17'], r['description17'])
        out.append(dict(r, x=math.log2(enr), y=-math.log10(max(p, 1e-300)),
                        p=p, cat=cat, color=color))
    return out


Y_BREAK = (50.0, 60.0)      # empty y band that is skipped by the axis break
POINT_S = 13.0              # marker area (pt^2) in the scan


def _volcano_panel(ax, rows, s=POINT_S, annotate=True, legend=True, yrange=None):
    """Draw the scan into one axes (called once per y segment when broken).

    yrange limits which points are emitted, so a broken axis does not write the
    full point set into both segments; the legend labels always use the full
    per-category counts.
    """
    def keep(r):
        return yrange is None or (yrange[0] <= r['y'] < yrange[1])
    for cat, color, _ in CATS + [('annotated (other)', '#BBBBBB', ''),
                                 ('unannotated', '#E6E6E6', '')]:
        sub = [r for r in rows if r['cat'] == cat]
        sel = [r for r in sub if keep(r)]
        if not sel:
            continue
        ax.scatter([r['x'] for r in sel], [r['y'] for r in sel], s=s, c=color,
                   alpha=0.6, lw=0,
                   label=f'{cat} ({len(sub)})' if legend else None)
    ax.axvline(0, color='#888888', lw=0.6, ls='--')
    ax.axhline(-math.log10(0.05), color='#888888', lw=0.6, ls=':')
    ax.axvline(-math.log2(1.3), color='#BBBBBB', lw=0.5, ls='-.')
    ax.axvline(math.log2(1.3), color='#BBBBBB', lw=0.5, ls='-.')
    if annotate:
        for r in rows:
            if r['og'] in HIGHLIGHT and keep(r):
                dx, dy = OFFSETS.get(r['og'], (3, 2))
                ax.annotate(HIGHLIGHT[r['og']], (r['x'], r['y']), xytext=(dx, dy),
                            textcoords='offset points',
                            ha='left' if dx > 0 else 'right', va='bottom',
                            fontsize=5.8, color='#7B2C29', zorder=20,
                            bbox=dict(fc='white', ec='none', alpha=0.6, pad=0.1),
                            arrowprops=dict(arrowstyle='-', lw=0.4, color='#999999',
                                            shrinkA=0, shrinkB=1))
    ax.tick_params(labelsize=6.4)
    ax.margins(x=0.03)
    return ax.get_legend_handles_labels()


def _break_marks(ax_lo, size_pt=2.4, gap_pt=1.6, color='#3C3C3C'):
    """Two 45-degree slashes on the spine at the axis break.

    The spine is interrupted between the slashes (white patch) so the marks
    read as a real break.  Sizes are given in points and converted to figure
    coordinates, so the marks are identical whatever the segment heights are.
    Only one pair is drawn: the two segments touch, and a pair on each side
    would overlap into a cross.
    """
    fig = ax_lo.figure
    fw = fig.get_size_inches()[0] * 72.0
    fh = fig.get_size_inches()[1] * 72.0
    pos = ax_lo.get_position()
    x0, y0 = pos.x0, pos.y1
    d = size_pt / fw                     # half width  of one slash
    h = size_pt / fh                     # half height of one slash
    g = gap_pt / fh                      # half gap between the two slashes
    ax_lo.add_patch(Rectangle((x0 - 2.2 * d, y0 - g - 1.2 * h), 4.4 * d,
                              2 * (g + 1.2 * h), transform=fig.transFigure,
                              clip_on=False, fc='white', ec='none', zorder=29))
    kw = dict(transform=fig.transFigure, color=color, clip_on=False, lw=0.9,
              solid_capstyle='butt', zorder=30)
    for k in (+1, -1):
        yc = y0 + k * g
        ax_lo.plot([x0 - d, x0 + d], [yc - h, yc + h], **kw)


def fig2c_broken(ax_lo, ax_hi, s=POINT_S, annotate=True, legend=True,
                 yticks_lo=(0, 10, 20, 30, 40), yticks_hi=(60, 80, 100, 120)):
    """Volcano with a broken y axis; returns (handles, labels)."""
    rows = load_points()
    mid = (Y_BREAK[0] + Y_BREAK[1]) / 2.0
    h, l = _volcano_panel(ax_lo, rows, s=s, annotate=annotate, legend=legend,
                          yrange=(-1e9, mid))
    _volcano_panel(ax_hi, rows, s=s, annotate=annotate, legend=False,
                   yrange=(mid, 1e9))
    ax_lo.set_ylim(0, Y_BREAK[0])
    ax_hi.set_ylim(Y_BREAK[1], 120)
    ax_hi.set_yticks(list(yticks_hi))
    ax_lo.set_yticks(list(yticks_lo))
    for t in list(ax_hi.get_yticklabels()) + list(ax_lo.get_yticklabels()):
        t.set_clip_on(False)      # keep the ticks that sit on the seam visible
    ax_hi.spines['bottom'].set_visible(False)
    ax_lo.spines['top'].set_visible(False)
    ax_hi.tick_params(bottom=False, labelbottom=False)
    _break_marks(ax_lo)
    return h, l


def fig2c(ax=None, figsize=(5.8, 4.3)):
    rows = load_points()
    own = ax is None
    if own:
        fig = plt.figure(figsize=(5.7, 4.1))
        outer = fig.add_gridspec(2, 2, width_ratios=[4.2, 0.95],
                                 height_ratios=[0.95, 4.4], wspace=0.07, hspace=0.07,
                                 left=0.105, right=0.985, top=0.94, bottom=0.155)
        ax_top = fig.add_subplot(outer[0, 0])
        vol = outer[1, 0].subgridspec(2, 1, height_ratios=[0.52, 2.33], hspace=0.0)
        ax_hi = fig.add_subplot(vol[0, 0])
        ax = fig.add_subplot(vol[1, 0], sharex=ax_hi)
        right = outer[1, 1].subgridspec(2, 1, height_ratios=[0.52, 2.33], hspace=0.0)
        ax_r_hi = fig.add_subplot(right[0, 0])
        ax_r_hi.axis('off')
        ax_right = fig.add_subplot(right[1, 0], sharey=ax)
        ax_corner = fig.add_subplot(outer[0, 1])
        ax_corner.axis('off')
        handles, labels = fig2c_broken(ax, ax_hi)
        ax.set_xlabel('log2 enrichment (racemose vs diploid Prunus)', fontsize=7.2)
        ax.set_ylabel('\u2212log10 p', fontsize=7.2)
    else:
        handles, labels = _volcano_panel(ax, rows)
        ax.set_ylim(0, 120)
    if own:
        # marginal distributions, split at x = 0 to mirror the two wings
        xs = [r['x'] for r in rows]
        ys = [r['y'] for r in rows]
        step = 0.1
        b0 = math.floor(min(xs) / step) * step
        b1 = math.ceil(max(xs) / step) * step
        bins_x = [b0 + i * step for i in range(int(round((b1 - b0) / step)) + 1)]
        edges_x = np.asarray(bins_x)
        cnt_l, _ = np.histogram([v for v in xs if v < 0], bins=bins_x)
        cnt_r, _ = np.histogram([v for v in xs if v >= 0], bins=bins_x)
        # bars start at the panel floor instead of y=0: on a log axis y=0 is -inf
        # and matplotlib would emit absurd coordinates (~1e4 pt) for the bar bottoms
        pos = np.concatenate([cnt_l[cnt_l > 0], cnt_r[cnt_r > 0]])
        lo, hi = pos.min(), pos.max()
        span = math.log10(hi) - math.log10(lo)
        base = 10 ** (math.log10(lo) - 0.05 * span)   # = autoscaled ylim bottom
        ax_top.set_ylim(base, 10 ** (math.log10(hi) + 0.05 * span))
        for cnt, col in ((cnt_l, '#9DB8D6'), (cnt_r, '#E0B0A8')):
            m = cnt > 0
            ax_top.bar(edges_x[:-1][m], cnt[m] - base, width=step, bottom=base,
                       align='edge', color=col, lw=0)
        ax_top.set_yscale('log')
        ax_top.set_yticks([])
        ax_top.spines['left'].set_visible(False)
        ax_top.axvline(0, color='#888888', lw=0.5, ls='--')
        ax_top.tick_params(bottom=False, labelbottom=False, which='both')
        plt.setp(ax_top.get_xticklabels(), visible=False)
        plt.setp(ax_top.get_xticklines(), visible=False)
        cnt, edges = np.histogram(ys, bins=72)
        # only bins inside the visible y range: the rest would be clipped away
        # but still carry coordinates far above the panel
        m = (cnt > 0) & (edges[:-1] < Y_BREAK[0])
        ax_right.barh(edges[:-1][m], cnt[m] - 0.7, height=edges[1] - edges[0],
                      left=0.7, align='edge', color='#CFCFCF', lw=0)
        ax_right.set_xscale('log')
        ax_right.set_xlim(0.7, max(cnt.max() * 2.0, 10))
        ax_right.set_xticks([])
        ax_right.tick_params(labelleft=False, left=False)
        ax_right.spines['bottom'].set_visible(False)

        n_all = len(rows)
        n_up = sum(1 for r in rows if r['x'] > 0 and float(r['q_cs']) < 0.05)
        n_dn = sum(1 for r in rows if r['x'] < 0 and float(r['q_cs_less']) < 0.05)
        n_cons = sum(1 for r in rows if float(r['q_cs']) < 0.05
                     and float(r['enr_cs']) > 1.3 and int(r['n_species_up']) >= 2
                     and r['te17'] != 'True')
        n_te = sum(1 for r in rows if r['cat'] == 'TE / transposon')
        ax_corner.text(0.02, 0.99,
                       f'{n_all:,} families tested\n'
                       f'{n_up:,} retained \u00b7 {n_dn:,} fractionated  (FDR < 0.05)\n'
                       f'{n_cons} consistent  (>1.3\u00d7, \u22652/3 spp., non-TE)\n'
                       f'{n_te} TE-flagged (orange)',
                       transform=ax_corner.transAxes, fontsize=4.9, va='top',
                       color='#555555', linespacing=1.5)
        ax.legend(handles, labels, fontsize=5.2, loc='upper center',
                  bbox_to_anchor=(0.5, -0.165), ncol=3, columnspacing=0.8,
                  handletextpad=0.25, labelspacing=0.35)
        ax_top.set_title('b  Genome-wide family scan (annotation from all 17 tracks)',
                         fontsize=8.4, loc='left')
        for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
            fig.savefig(FIG / f'fig2c_scan_volcano_all17.{ext}', bbox_inches='tight', **kw)
        plt.close(fig)
        counts = Counter(r['cat'] for r in rows)
        print('fig2c saved; plotted', len(rows), dict(counts))
    else:
        ax.legend(handles, labels, fontsize=4.0, loc='upper left', ncol=2,
                  columnspacing=0.6, handletextpad=0.2, labelspacing=0.22)
    return ax


if __name__ == '__main__':
    fig2c()
