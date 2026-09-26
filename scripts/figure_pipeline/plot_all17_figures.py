#!/usr/bin/env python3
"""plot_all17_figures.py - all-17-annotation versions of the GO figures and the
integrated summary figure.

  fig4_all17 : GO enrichment of the racemose-retained set (all-17 annotation)
  fig7_all17 : overlap of the three clade sets + GO profiles (all-17 annotation,
               includes the shared corymbose & solitary set)
  fig0_all17 : integrated six-panel summary (a-f) with the all-17 volcano
               (panel c), the all-17 GO panel (panel e) and the set-overlap /
               GO-profile panel (panel f); other panels reuse
               plot_figures_v2.py unchanged.

Outputs: fig4_go_enrichment_all17.*, fig7_shared_specific_all17.*,
         fig3_top_family_heatmap.*, fig0_integrated.*
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
from matplotlib.patches import Circle

import plot_figures_v2 as pv
from plot_fig2c_scan_all17 import fig2c, classify
from plot_fig4b_go_clades import go_bubble_panel, size_legend_handles
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
RES = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'


def save(fig, name):
    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        fig.savefig(FIG / f'{name}.{ext}', bbox_inches='tight', **kw)
    plt.close(fig)
    print('saved', name)


# ---------------------------------------------------------------- fig4 (all17)
def fig4_all17(ax=None):
    rows = [r for r in csv.DictReader((RES / 'go_enrichment_racemose_all17.tsv').open(),
                                      delimiter='\t')]
    rows = [r for r in rows if float(r['q']) < 0.25][:8][::-1]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(5.0, 3.4))
    go_bubble_panel(ax, rows, '#C0504D', label_len=42)
    if own:
        handles = size_legend_handles('#8C8C8C')
        fig.legend(handles=handles, title='families with term', fontsize=5.4,
                   title_fontsize=5.6, loc='lower right', bbox_to_anchor=(1.0, -0.01),
                   ncol=3, columnspacing=1.4, handletextpad=0.6)
        fig.text(0.01, 0.005,
                 'bubble area \u221d \u221a(families with term); label = '
                 'families/background; dotted line: FDR = 0.05',
                 fontsize=5.2, color='#555555')
        ax.set_title('d  GO enrichment (racemose-retained; all-17 annotation)',
                     fontsize=8.4, loc='left')
        save(fig, 'fig4_go_enrichment_all17')
    return ax


# ---------------------------------------------------------------- fig7 (all17)
def fig7_all17(axs=None, sub_titles=True):
    counts = {}
    with (RES / 'shared_specific_counts_all17.tsv').open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['category']] = int(r['n_families'])
    rows = list(csv.DictReader((RES / 'shared_specific_go_all17.tsv').open(), delimiter='\t'))
    own = axs is None
    if own:
        fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.2),
                                 gridspec_kw={'width_ratios': [0.95, 1.35], 'wspace': 0.35})
    else:
        axes = list(axs)
    axa = axes[0]
    axa.set_xlim(-1.55, 1.55)
    axa.set_ylim(-1.45, 1.25)
    axa.set_aspect('equal')          # keep the three circles round
    axa.axis('off')
    for (cx, cy), col in [((-0.42, 0.28), '#C0504D'), ((0.42, 0.28), '#4F81BD'),
                          ((0.0, -0.42), '#59A14F')]:
        axa.add_patch(Circle((cx, cy), 0.78, fc=col, alpha=0.20, ec=col, lw=1.0))
    put = lambda x, y, n, fs=7.2: axa.text(x, y, str(n), ha='center', va='center', fontsize=fs)
    put(-0.92, 0.32, counts['R_only'])
    put(0.92, 0.32, counts['C_only'])
    put(0.0, -0.92, counts['S_only'])
    put(0.0, 0.58, counts['R_C'], 6.4)
    put(-0.60, -0.50, counts['R_S'], 6.4)
    put(0.60, -0.50, counts['C_S'], 6.4)
    put(0.0, -0.16, counts['R_C_S'], 6.4)
    axa.text(-1.05, 1.05, f'racemose\n(n={counts["R_only"] + counts["R_C"] + counts["R_S"] + counts["R_C_S"]:,})', color='#C0504D', fontsize=6.4, ha='center')
    axa.text(1.05, 1.05, f'corymbose\n(n={counts["C_only"] + counts["R_C"] + counts["C_S"] + counts["R_C_S"]:,})', color='#4F81BD', fontsize=6.4, ha='center')
    axa.text(0.0, -1.32, f'solitary (n={counts["S_only"] + counts["R_S"] + counts["C_S"] + counts["R_C_S"]:,})', color='#59A14F', fontsize=6.4, ha='center')
    if own and sub_titles:
        axa.set_title('a  Overlap of the enrichment sets', fontsize=8.0, loc='left')

    axb = axes[1]
    CAT_COL = {'R_only': '#C0504D', 'C_only': '#4F81BD', 'S_only': '#59A14F',
               'C_S': '#9A4D8E'}
    picks = [('R_only', 4), ('C_only', 4), ('S_only', 4), ('C_S', 3)]
    sel = []
    for g, n in picks:
        rs = sorted([r for r in rows if r['group'] == g], key=lambda r: float(r['q']))[:n]
        sel += rs[::-1]
    ys = np.arange(len(sel))
    for y, r in zip(ys, sel):
        col = CAT_COL[r['group']]
        x = -math.log10(max(float(r['q']), 1e-12))
        axb.plot([0, x], [y, y], color=col, lw=0.8, alpha=0.55)
        axb.scatter([x], [y], s=10 + 4 * int(r['k']), color=col, zorder=3)
    axb.set_yticks(ys, [r['go_description'][:38] for r in sel], fontsize=5.2)
    for tick, r in zip(axb.get_yticklabels(), sel):
        tick.set_color(CAT_COL[r['group']])
    axb.set_xlabel('\u2212log10 q (BH)', fontsize=7.0)
    axb.tick_params(axis='x', labelsize=6.4)
    axb.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    axb.set_xlim(0, max(-math.log10(max(float(r['q']), 1e-12)) for r in sel) * 1.12)
    handles = [Line2D([], [], marker='o', ls='', color=c, label=l, ms=4) for l, c in
               [('racemose-only', '#C0504D'), ('corymbose-only', '#4F81BD'),
                ('solitary-only', '#59A14F'), ('shared cor & sol', '#9A4D8E')]]
    axb.legend(handles=handles, fontsize=5.4, loc='lower right')
    if own:
        if sub_titles:
            axb.set_title('b  GO profiles of the clade sets (all-17 annotation)',
                          fontsize=8.0, loc='left')
        save(fig, 'fig7_shared_specific_all17')
    return axes


SPECIES_COLORS = {
    'pzi_A': '#B64342', 'pzi_B': '#B64342', 'psp_A': '#E28E2C', 'psp_B': '#E28E2C',
    'ppa_A': '#9A4D8E', 'ppa_B': '#9A4D8E',
    'pcam': '#5B8FD6', 'pav': '#5B8FD6', 'pspe': '#5B8FD6', 'pmal': '#5B8FD6',
    'par': '#7BAA5B', 'pmu': '#7BAA5B', 'pdu': '#7BAA5B', 'psal': '#7BAA5B',
    'ppe': '#7BAA5B', 'gtr': '#8F8F8F', 'rch': '#8F8F8F',
}


def _top_strip(ax, tracks, y0=-1.02, h=0.42):
    for j, t in enumerate(tracks):
        ax.add_patch(Rectangle((j - 0.5, y0), 1.0, h, color=SPECIES_COLORS[t],
                               clip_on=False, lw=0))


# ---------------------------------------------------------------- fig3 (all17)
def fig3_all17(ax=None):
    rows = list(csv.DictReader((FIG / 'data/clade_family_scan_v3.tsv').open(), delimiter='\t'))
    counts = {}
    with (B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups/'
              'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['Orthogroup']] = r
    tracks = ['pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B',
              'pcam', 'pav', 'pspe', 'pmal', 'par', 'pmu', 'pdu', 'psal', 'ppe', 'gtr', 'rch']
    glabels = ['PZI-A', 'PZI-B', 'PSP-A', 'PSP-B', 'PPA-A', 'PPA-B',
               'PCAM', 'PAV', 'PSPE', 'PMAL', 'PAR', 'PMU', 'PDU', 'PSAL', 'PPE', 'GTR', 'RCH']
    sel = [r for r in rows
           if (r['symbol17'] or '').strip() and float(r['rac']) >= 7
           and float(r['p_cs']) < 0.05 and float(r['enr_cs']) > 1.3
           and int(r['n_species_up']) >= 2 and r['te17'] != 'True']
    sel.sort(key=lambda r: -float(r['enr_cs']))
    top = sel[:12]
    M = np.zeros((len(top), len(tracks)))
    cats = []
    for i, r in enumerate(top):
        cr = counts[r['og']]
        for j, t in enumerate(tracks):
            M[i, j] = int(cr[t])
        cats.append(classify(r['symbol17'], r['description17']))
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
    im = ax.imshow(np.log1p(M), cmap='Reds', aspect='auto', vmin=0, vmax=np.log1p(M).max())
    ax.set_xticks(range(len(tracks)), glabels, rotation=90, fontsize=5.8)
    ax.set_yticks(range(len(top)),
                  [f"{r['symbol17'].rstrip(';')} ({r['og'][-5:]}) "
                   + ('***' if int(r['n_species_up']) >= 3 else '**') for r in top],
                  fontsize=6.0)
    _top_strip(ax, tracks)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = int(M[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center', fontsize=5.0,
                        color='white' if v > 6 else 'black')
    ax.axvline(5.5, color='black', lw=0.9)
    ax.axvline(9.5, color='#666666', lw=0.6, ls='--')
    ax.axvline(14.5, color='#666666', lw=0.6, ls='--')
    ax.text(2.5, -1.15, 'racemose', fontsize=6.4, ha='center', color='#C0504D')
    ax.text(7.5, -1.15, 'corymbose clade', fontsize=6.4, ha='center', color='#555555')
    ax.text(12.0, -1.15, 'solitary clade', fontsize=6.4, ha='center', color='#555555')
    ax.text(15.5, -1.15, 'outgroup', fontsize=6.4, ha='center', color='#555555')
    code = {'TE / transposon': 'TE', 'floral / meristem / hormone': 'floral',
            'defense / stress': 'defense', 'RNA / splicing': 'RNA',
            'transport / vesicle': 'transport', 'meiosis / DNA / chromatin': 'meiosis',
            'metabolism / other': 'metab', 'annotated (other)': 'other',
            'unannotated': '-'}
    for i, (cat, color) in enumerate(cats):
        ax.add_patch(Rectangle((len(tracks) - 0.35, i - 0.5), 0.5, 1.0, color=color,
                               clip_on=False))
    ax.set_xlim(-0.5, len(tracks) + 0.4)
    used = []
    for cat, color in cats:
        if cat not in [u[0] for u in used]:
            used.append((cat, color))
    handles = [Line2D([], [], marker='s', ls='', color=c, label=code[k], ms=5)
               for k, c in used]
    ax.legend(handles=handles, fontsize=5.4, loc='lower left',
              bbox_to_anchor=(0.62, 0.02), ncol=2, columnspacing=0.8,
              handletextpad=0.3)
    if own:
        ax.set_title('c  Consistent racemose-retained families '
                     '(copy numbers; right = category; all-17 annotation)',
                     fontsize=8.4, loc='left', pad=16)
        save(fig, 'fig3_top_family_heatmap')
    return ax


# ---------------------------------------------------------------- fig3 detailed
def fig3_detailed():
    rows = list(csv.DictReader((FIG / 'data/clade_family_scan_v3.tsv').open(), delimiter='\t'))
    counts = {}
    with (B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups/'
              'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['Orthogroup']] = r
    tracks = ['pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B',
              'pcam', 'pav', 'pspe', 'pmal', 'par', 'pmu', 'pdu', 'psal', 'ppe', 'gtr', 'rch']
    glabels = ['PZI-A', 'PZI-B', 'PSP-A', 'PSP-B', 'PPA-A', 'PPA-B',
               'PCAM', 'PAV', 'PSPE', 'PMAL', 'PAR', 'PMU', 'PDU', 'PSAL', 'PPE', 'GTR', 'RCH']
    sel = [r for r in rows
           if (r['symbol17'] or '').strip() and float(r['rac']) >= 7
           and float(r['p_cs']) < 0.05 and float(r['enr_cs']) > 1.3
           and int(r['n_species_up']) >= 2 and r['te17'] != 'True']
    sel.sort(key=lambda r: -float(r['enr_cs']))
    top = sel[:12]
    n = len(top)
    M = np.zeros((n, len(tracks)))
    cats = []
    for i, r in enumerate(top):
        cr = counts[r['og']]
        for j, t in enumerate(tracks):
            M[i, j] = int(cr[t])
        cats.append(classify(r['symbol17'], r['description17']))

    fig = plt.figure(figsize=(10.2, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.32, 0.28], wspace=0.45,
                          left=0.145, right=0.965, top=0.77, bottom=0.30)
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_enr = fig.add_subplot(gs[0, 1])
    ax_sig = fig.add_subplot(gs[0, 2])

    # heatmap
    im = ax_hm.imshow(np.log1p(M), cmap='Reds', aspect='auto', vmin=0, vmax=np.log1p(M).max())
    ax_hm.set_xticks(range(len(tracks)), glabels, rotation=90, fontsize=5.8)
    ax_hm.set_yticks(range(n),
                     [f"{r['symbol17'].rstrip(';')} ({r['og'][-5:]}) "
                      + ('***' if int(r['n_species_up']) >= 3 else '**') for r in top],
                     fontsize=6.2)
    ax_hm.set_xticks(np.arange(-0.5, len(tracks), 1), minor=True)
    ax_hm.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax_hm.grid(which='minor', color='white', lw=0.6)
    ax_hm.tick_params(which='minor', length=0)
    for i in range(n):
        for j in range(len(tracks)):
            v = int(M[i, j])
            if v > 0:
                ax_hm.text(j, i, str(v), ha='center', va='center', fontsize=5.0,
                           color='white' if v > 7 else 'black')
    ax_hm.axvline(5.5, color='black', lw=1.0)
    ax_hm.axvline(9.5, color='#666666', lw=0.6, ls='--')
    ax_hm.axvline(14.5, color='#666666', lw=0.6, ls='--')
    _top_strip(ax_hm, tracks, y0=-1.06, h=0.40)
    ax_hm.text(2.5, -1.62, 'racemose', fontsize=6.4, ha='center', color='#C0504D')
    ax_hm.text(7.5, -1.62, 'corymbose clade', fontsize=6.4, ha='center', color='#555555')
    ax_hm.text(12.0, -1.62, 'solitary clade', fontsize=6.4, ha='center', color='#555555')
    ax_hm.text(15.5, -1.62, 'outgroup', fontsize=6.4, ha='center', color='#555555')

    # A/B copy split drawn inside the heatmap axes at negative x
    ab_tot = [float(r['racA']) + float(r['racB']) for r in top]
    ab_scale = 1.30 / max(ab_tot)
    for i, r in enumerate(top):
        a, b = float(r['racA']), float(r['racB'])
        ax_hm.barh(i, a * ab_scale, left=-2.30, color='#B64342', height=0.68, zorder=3)
        ax_hm.barh(i, b * ab_scale, left=-2.30 + a * ab_scale, color='#4F81BD',
                   height=0.68, zorder=3)
        ax_hm.text(-2.30 + (a + b) * ab_scale + 0.05, i, f'{int(a + b)}', va='center',
                   fontsize=5.0, color='#666666', zorder=3)
    ax_hm.text(-1.42, -1.62, 'A / B copies', fontsize=6.0, ha='center', color='#555555')
    ax_hm.set_xlim(-2.55, len(tracks) - 0.5)
    for i, (cat, color) in enumerate(cats):
        ax_hm.add_patch(Rectangle((len(tracks) - 0.42, i - 0.5), 0.32, 1.0,
                                  color=color, clip_on=False, lw=0))

    # right: enrichment and significance
    enr = [float(r['enr_cs']) for r in top]
    le = [math.log2(v) for v in enr]
    ax_enr.barh(range(n), le, color='#C0504D', height=0.68, alpha=0.9)
    for i, v in enumerate(enr):
        ax_enr.text(math.log2(v) + max(le) * 0.04, i, f'{v:.1f}\u00d7', va='center',
                    fontsize=5.0, color='#7B2C29')
    ax_enr.set_ylim(n - 0.5, -0.5)
    ax_enr.set_xlim(0, max(le) * 1.5)
    ax_enr.set_xticks([]); ax_enr.set_yticks([])
    ax_enr.set_title('enrichment', fontsize=6.0, pad=6, color='#555555')
    for s in ('top', 'right', 'bottom'):
        ax_enr.spines[s].set_visible(False)

    sig = [-math.log10(max(float(r['p_cs']), 1e-300)) for r in top]
    cap = 35.0
    ax_sig.barh(range(n), [min(s, cap) for s in sig], color='#4F81BD', height=0.68, alpha=0.9)
    for i, r, s in zip(range(n), top, sig):
        ax_sig.text(min(s, cap) + max(sig) * 0.03, i, f"{float(r['p_cs']):.0e}",
                    va='center', fontsize=5.0, color='#2C4E7B')
    ax_sig.set_ylim(n - 0.5, -0.5)
    ax_sig.set_xlim(0, max(sig) * 1.42)
    ax_sig.set_xticks([]); ax_sig.set_yticks([])
    ax_sig.set_title('\u2212log10 p', fontsize=6.0, pad=6, color='#555555')
    for s in ('top', 'right', 'bottom'):
        ax_sig.spines[s].set_visible(False)

    pos = ax_hm.get_position()
    cax = fig.add_axes([pos.x0 + 0.08 * pos.width, pos.y0 - 0.155,
                        0.62 * pos.width, 0.022])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_ticks(np.log1p([0, 1, 2, 5, 10, 25, 50, 100]))
    cb.set_ticklabels(['0', '1', '2', '5', '10', '25', '50', '100'])
    cb.ax.tick_params(labelsize=5.6)
    cb.set_label('gene copies per track', fontsize=6.0)

    code = {'TE / transposon': 'TE', 'floral / meristem / hormone': 'floral',
            'defense / stress': 'defense', 'RNA / splicing': 'RNA',
            'transport / vesicle': 'transport', 'meiosis / DNA / chromatin': 'meiosis',
            'metabolism / other': 'metab', 'annotated (other)': 'other',
            'unannotated': '-'}
    used = []
    for cat, color in cats:
        if cat not in [u[0] for u in used]:
            used.append((cat, color))
    handles = [Line2D([], [], marker='s', ls='', color=c, label=code[k], ms=5)
               for k, c in used]
    fig.legend(handles=handles, fontsize=5.4, loc='lower left',
               bbox_to_anchor=(0.145, 0.042), ncol=4, columnspacing=0.9,
               handletextpad=0.3, frameon=False)
    fig.text(0.145, 0.008,
             'right panels: enrichment (log2) and significance (\u2212log10 p); stars: '
             'consistent in 2/3 (**) or 3/3 (***) racemose species; left bar: A (red) / '
             'B (blue) copies; top strip: clade colours; colourbar: copies per track '
             '(log1p scale)',
             fontsize=5.2, color='#555555')
    fig.suptitle('c  Consistent racemose-retained families (copies, enrichment, '
                 'significance; all-17 annotation)', fontsize=8.6, x=0.145,
                 ha='left', y=0.955)
    save(fig, 'fig3_top_family_heatmap')


# ---------------------------------------------------------------- fig0 (integrated, all17)
def fig0_all17():
    fig = plt.figure(figsize=(11.5, 9.9))
    gs = fig.add_gridspec(3, 3, hspace=0.46, wspace=0.40,
                          left=0.06, right=0.985, top=0.945, bottom=0.055,
                          width_ratios=[1.0, 0.98, 0.95],
                          height_ratios=[1.0, 1.05, 0.66])

    axa = fig.add_subplot(gs[0, 0])
    root, H, ids = pv.load_tree()
    pv.draw_chronogram(axa, root, H, tags={'<11>': ('WGD', '#B03A2E')},
                       xpad=46, fs_label=5.4, show_axis=False)
    for node in pv.iter_nodes(root):
        nm = node['name']
        if nm in pv.INFL_TYPE:
            pv.icon_at(axa, H + 35, node['y'], 0.22, pv.INFL_TYPE[nm])
    axa.set_title('a  Framework: dated tree, inflorescence types and the WGD',
                  fontsize=8.4, loc='left')

    axb = fig.add_subplot(gs[0, 1])
    pv.fig1(axb, compact=True)
    axb.set_title('b  CAFE5 family-size change along the tree', fontsize=8.4, loc='left')

    axc = fig.add_subplot(gs[0, 2])
    fig2c(axc)
    axc.set_title('c  Genome-wide family scan (all 17 tracks annotated)',
                  fontsize=8.4, loc='left')

    axd = fig.add_subplot(gs[1, 0:2])
    fig3_all17(axd)
    axd.set_title('d  Consistent racemose-retained families (right: functional category)',
                  fontsize=8.4, loc='left', pad=16)

    axe = fig.add_subplot(gs[1, 2])
    fig4_all17(axe)
    axe.set_title('e  GO enrichment of retained families (all-17 annotation)',
                  fontsize=8.4, loc='left')

    gsf = gs[2, 0:3].subgridspec(1, 2, width_ratios=[0.95, 1.35], wspace=0.30)
    axf1 = fig.add_subplot(gsf[0, 0])
    axf2 = fig.add_subplot(gsf[0, 1])
    fig7_all17(axs=(axf1, axf2), sub_titles=False)
    axf1.set_title('f  Shared vs clade-specific expansions: set overlap',
                   fontsize=8.4, loc='left')
    axf2.set_title('GO profiles of the clade sets', fontsize=8.4, loc='left')

    save(fig, 'fig0_integrated')


if __name__ == '__main__':
    fig3_detailed()
    fig4_all17()
    fig7_all17()
    fig0_all17()
