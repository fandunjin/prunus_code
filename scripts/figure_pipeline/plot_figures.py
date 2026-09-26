#!/usr/bin/env python3
"""plot_figures.py - figure set for the racemose-specific family scan.

fig1  CAFE5 branch landscape (expansion vs contraction per branch)
fig2  genome-wide scan volcano (racemose vs diploid Prunus)
fig3  top consistent racemose-expanded families across the 17 tracks
fig4  GO enrichment of consistent racemose-expanded families
fig5  synthesis model of inflorescence evolution

All figures are written as png/pdf/svg to 0903/figures/inflorescence_evolution_review/
"""
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from scipy.stats import binomtest

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

B = Path('${PROJECT_ROOT}')
SCAN = B / 'analysis/racemose_inflorescence_scan_20260922'
FIG = B / '0903/figures/inflorescence_evolution_review'
FIG.mkdir(parents=True, exist_ok=True)
PREV = Path('${HOME}/AppData/Local/Temp/opencode/fig1e_preview')

TOT = {'RAC': 186830, 'RAC_A': 93415, 'RAC_B': 93415, 'PZI': 62277, 'PSP': 93853,
       'PPA': 56318, 'COR': 131314, 'SOL': 150034, 'OUT': 56748}


def save(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(FIG / f'{name}.{ext}', bbox_inches='tight')
    fig.savefig(PREV / f'{name}_preview.png', dpi=150, bbox_inches='tight')
    print('saved', name)


# ------------------------------------------------------------------ fig1
def fig1():
    rows = []
    with (B / '0903/manuscript_pipeline/stage2_write/phase4_draft/figure5/output/source_data/'
              'cafe5_clade_branch_summary.tsv').open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['branch'] in ('<8>', '<9>', '<10>', '<11>', '<12>', '<13>'):
                continue
            rows.append(r)
    labelmap = {'<11>': 'racemose stem (WGD)', '<12>': 'psp–pzi stem',
                '<10>': 'stem (racemose+pcam/ppe)', '<13>': 'pcam–ppe stem',
                '<9>': 'gtr–Prunus stem'}
    order = ['<11>', '<12>', '<10>', '<13>', '<9>', 'PZI', 'PSP', 'PPA', 'PCAM', 'PPE',
             'GTR', 'RCH']
    by = {r['branch']: r for r in rows}
    labels, ups, downs = [], [], []
    for b in order:
        r = by.get(b)
        if not r:
            continue
        labels.append(labelmap.get(b, b))
        ups.append(int(r['families_changed_up']))
        downs.append(int(r['families_changed_down']))
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    ax.barh(y, ups, color='#C0504D', height=0.66, label='expanded')
    ax.barh(y, [-d for d in downs], color='#4F81BD', height=0.66, label='contracted')
    ax.axvline(0, color='black', lw=0.6)
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_xlabel('CAFE5 families with significant change', fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6.5, loc='lower right')
    for yy, u, d in zip(y, ups, downs):
        ax.text(u + 120, yy, f'{u:,}', va='center', fontsize=5.8, color='#7B2C29')
        ax.text(-d - 120, yy, f'{d:,}', va='center', ha='right', fontsize=5.8, color='#2C4E7B')
    ax.set_xlim(-3400, 10500)
    ax.set_title('a  Family-size change per branch', fontsize=8.4, loc='left')
    save(fig, 'fig1_cafe5_branch_landscape')


# ------------------------------------------------------------------ fig2
def fig2():
    rows = list(csv.DictReader((SCAN / 'results_v2/clade_family_scan.tsv').open(), delimiter='\t'))
    xs, ys, colors, labels = [], [], [], []
    for r in rows:
        if float(r['rac']) + float(r['cor']) + float(r['sol']) + float(r['out']) == 0:
            continue
        enr = float(r['enr_cs']) if r['enr_cs'] not in ('', 'nan') else float('nan')
        p_up = float(r['p_cs'])
        base = float(r['cor']) + float(r['sol'])
        p_dn = binomtest(int(r['rac']), TOT['RAC'],
                         min(1.0, base / (TOT['COR'] + TOT['SOL'])),
                         alternative='less').pvalue if base > 0 else 1.0
        if not math.isfinite(enr) or enr <= 0:
            continue
        x = math.log2(enr)
        p = min(p_up, p_dn)
        y = -math.log10(max(p, 1e-300))
        xs.append(x)
        ys.append(y)
        te = r['te'] == 'True'
        sig = p < 0.05 and abs(x) > math.log2(1.3)
        colors.append('#D9D9D9' if not sig else ('#E8A33D' if te else '#C0504D'))
        labels.append(r)
    fig, ax = plt.subplots(figsize=(4.2, 3.1))
    ax.scatter(xs, ys, s=3.5, c=colors, alpha=0.65, lw=0)
    ax.axvline(0, color='#888888', lw=0.6, ls='--')
    ax.axhline(-math.log10(0.05), color='#888888', lw=0.6, ls=':')
    HIGHLIGHT = {'OG0000961': 'BLH1', 'OG0005156': 'ASY1', 'OG0005416': 'PDR2',
                 'OG0000864': 'CPK1', 'OG0000130': 'LNK1', 'OG0000222': 'RAPTOR2',
                 'OG0001188': 'Susy', 'OG0000664': 'enolase', 'OG0000348': 'MCM7'}
    for r, x, y in zip(labels, xs, ys):
        if r['og'] in HIGHLIGHT:
            ax.annotate(HIGHLIGHT[r['og']], (x, y), xytext=(2, 2), textcoords='offset points',
                        fontsize=6.0, color='#7B2C29')
    ax.set_xlabel('log2 enrichment (racemose vs diploid Prunus)', fontsize=7.5)
    ax.set_ylabel('\u2212log10 p', fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.text(-3.6, 4.4, 'fractionated\nin racemose', fontsize=6, color='#4F81BD', ha='left')
    ax.text(2.0, 4.4, 'retained /\nexpanded', fontsize=6, color='#C0504D')
    ax.set_title('b  Genome-wide family scan', fontsize=8.4, loc='left')
    save(fig, 'fig2_scan_volcano')


# ------------------------------------------------------------------ fig3
def fig3():
    exp = list(csv.DictReader((SCAN / 'results_v2/racemose_expanded_consistent.tsv').open(), delimiter='\t'))
    counts, totals = {}, {}
    with (B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups/'
              'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['Orthogroup']] = r
    tracks = ['pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B',
              'pcam', 'pav', 'pspe', 'pmal', 'par', 'pmu', 'pdu', 'psal', 'ppe', 'gtr', 'rch']
    glabels = ['PZI-A', 'PZI-B', 'PSP-A', 'PSP-B', 'PPA-A', 'PPA-B',
               'PCAM', 'PAV', 'PSPE', 'PMAL', 'PAR', 'PMU', 'PDU', 'PSAL', 'PPE', 'GTR', 'RCH']
    top = [r for r in exp if r['symbol'] and float(r['rac']) >= 7][:22]
    M = np.zeros((len(top), len(tracks)))
    for i, r in enumerate(top):
        cr = counts[r['og']]
        for j, t in enumerate(tracks):
            M[i, j] = int(cr[t])
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    im = ax.imshow(np.log1p(M), cmap='Reds', aspect='auto', vmin=0, vmax=np.log1p(M).max())
    ax.set_xticks(range(len(tracks)), glabels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(top)),
                  [f"{r['symbol']}  ({r['og'][-5:]})" for r in top], fontsize=6.2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = int(M[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center', fontsize=5.2,
                        color='white' if v > 6 else 'black')
    ax.axvline(5.5, color='black', lw=0.9)
    ax.axvline(9.5, color='#666666', lw=0.6, ls='--')
    ax.axvline(14.5, color='#666666', lw=0.6, ls='--')
    ax.text(2.5, -0.75, 'racemose', fontsize=6.6, ha='center', color='#C0504D')
    ax.text(7.5, -0.75, 'corymbose clade', fontsize=6.6, ha='center', color='#555555')
    ax.text(12.0, -0.75, 'solitary clade', fontsize=6.6, ha='center', color='#555555')
    ax.text(15.5, -0.75, 'outgroup', fontsize=6.6, ha='center', color='#555555')
    ax.set_title('c  Consistent racemose-expanded families (copy numbers)', fontsize=8.4, loc='left', pad=16)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.ax.tick_params(labelsize=6)
    cb.set_label('log1p copies', fontsize=6.5)
    save(fig, 'fig3_top_family_heatmap')


# ------------------------------------------------------------------ fig4
def fig4():
    rows = [r for r in csv.DictReader((SCAN / 'results_v2/go_enrichment_consistent.tsv').open(),
                                      delimiter='\t')]
    rows = [r for r in rows if float(r['q']) < 0.25][:14]
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    y = np.arange(len(rows))
    sizes = [int(r['expanded']) for r in rows]
    ax.scatter([-math.log10(max(float(r['q']), 1e-12)) for r in rows], y,
               s=[18 + 8 * s for s in sizes], color='#C0504D', alpha=0.8, lw=0)
    for yy, r in zip(y, rows):
        ax.text(-math.log10(max(float(r['q']), 1e-12)) + 0.06, yy,
                f"{r['expanded']}/{r['bg']}", fontsize=5.8, va='center', color='#555555')
    ax.set_yticks(y, [r['go_description'][:42] for r in rows], fontsize=6.2)
    ax.set_xlabel('\u2212log10 FDR (BH)', fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    ax.set_title('d  GO enrichment (consistent racemose-expanded)', fontsize=8.4, loc='left')
    save(fig, 'fig4_go_enrichment')


# ------------------------------------------------------------------ fig5
def fig5():
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.axis('off')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.4)

    def box(x, y, w, h, text, fc, ec='#555555', fs=6.4):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.07',
                                    fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs)

    def arrow(x1, y1, x2, y2, color='#555555'):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                     mutation_scale=8, color=color, lw=0.9))

    box(0.15, 2.9, 2.0, 0.95, 'diploid Prunus\nancestor (n = 8)\nracemose flowers',
        '#F2F2F2')
    arrow(2.15, 3.35, 2.65, 3.35)
    box(2.65, 2.9, 1.9, 0.95, 'shared WGD\n(allopolyploidy)\n~55–42 Ma', '#FBE3D6')

    box(5.15, 4.15, 3.0, 1.75,
        'racemose lineage\n(PZI, PSP, PPA)\n\n'
        '9,727 families expanded\n'
        '466 retained above the\nWGD expectation\n'
        'SVP/DAM, SPL9/15, YABBY1/2/3,\nBLH1 (KNOX/BELL)', '#FBE3D6')
    arrow(4.55, 3.65, 5.15, 4.6)
    box(5.15, 2.05, 3.0, 1.25, 'corymbose clade\n(PCAM, PAV, PSPE, PMAL)\n\n'
        'contraction-dominated\n(2,244 families in PCAM)', '#E8EEF7')
    arrow(4.55, 3.3, 5.15, 2.9)
    box(5.15, 0.25, 3.0, 1.25, 'solitary clade\n(PPE, PSAL, PAR, PMU, PDU)\n\n'
        'contraction-dominated\n(1,489 families in PPE)', '#E8EEF7')
    arrow(4.55, 3.0, 5.15, 1.1)

    box(8.55, 4.25, 3.25, 1.55,
        'ancestral racemose\nprogramme retained\n\ncanalized by duplicated\nregulatory networks', '#DFF0D8')
    arrow(8.15, 5.0, 8.55, 5.0)
    box(8.55, 1.35, 3.25, 1.5,
        'derived inflorescence\narchitectures\n\ncopy-number reduction and\nregulatory divergence in\ndiploid descendants',
        '#E8EEF7')
    arrow(8.15, 2.65, 8.55, 2.2)

    box(5.15, 5.95, 6.65, 0.75,
        'non-Prunus outgroups (GTR, RCH): contraction-dominated\nno SVP/DAM expansion',
        '#FFFFFF', ec='#BBBBBB', fs=6.2)
    box(0.15, 5.95, 4.7, 0.75,
        'interpretation: ancestral racemose state retained\nin the polyploids; corymbose and solitary\nprogrammes diverged after the split',
        '#FFFFFF', ec='#BBBBBB', fs=5.8)
    ax.set_title('e  Working model: WGD-triggered retention of the ancestral racemose programme',
                 fontsize=8.4, loc='left')
    save(fig, 'fig5_evolution_model')


if __name__ == '__main__':
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
