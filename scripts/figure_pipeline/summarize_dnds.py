#!/usr/bin/env python3
"""summarize_dnds.py - integrate branch dN/dS results (candidates, controls,
genome-wide background A and B) and draw fig6.

Inputs : results/branch_dnds_results.tsv, results/branch_dnds_background.tsv,
         results/branch_dnds_background_B.tsv
Outputs: results/branch_dnds_integrated.tsv
         0903/figures/inflorescence_evolution_review/fig6_dnds_branch.{png,pdf,svg}
"""
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 7,
    'axes.spines.right': False, 'axes.spines.top': False, 'axes.linewidth': 0.8,
    'legend.frameon': False,
})

B = Path('${PROJECT_ROOT}')
D = B / 'analysis/dnds_branch_20260922'
RES = D / 'results'
FIG = B / '0903/figures/inflorescence_evolution_review'

CAND = {'SVP/DAM', 'KNAT1/2/6', 'STM', 'YABBY1/2/3', 'SPL9/15'}

rows = []
with (RES / 'branch_dnds_results.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        r['group'] = 'candidate' if r['subfamily'] in CAND else 'control'
        rows.append(r)
with (RES / 'branch_dnds_background.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        r['group'] = 'background_A'
        r['subfamily'] = ''
        r['family'] = 'single-copy OG'
        rows.append(r)
with (RES / 'branch_dnds_background_B.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        r['group'] = 'background_B'
        r['subfamily'] = ''
        r['family'] = 'single-copy OG'
        rows.append(r)

with (RES / 'branch_dnds_integrated.tsv').open('w', newline='') as fh:
    cols = ['group', 'family', 'subfamily', 'set', 'og', 'n_codons',
            'omega_bg', 'omega_rac', 'LRT', 'p']
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    for r in rows:
        w.writerow(r)

groups = defaultdict(list)
for r in rows:
    groups[r['group']].append(r)
print('group  n  p<0.05  mean w_rac  mean w_bg  ratio')
for g in ['candidate', 'control', 'background_A', 'background_B']:
    rs = groups[g]
    ps = [float(r['p']) for r in rs]
    wr = [float(r['omega_rac']) for r in rs if r['omega_rac'] not in ('', 'None')]
    wb = [float(r['omega_bg']) for r in rs if r['omega_bg'] not in ('', 'None')]
    sig = sum(1 for p in ps if p < 0.05)
    print(f'{g:14s} {len(rs):3d}  {sig:3d} ({sig/len(rs):.0%})  {np.mean(wr):.3f}  '
          f'{np.mean(wb):.3f}  {np.mean(wr)/np.mean(wb):.2f}')

# ------------------------------------------------------------------ figure
colors = {'candidate': '#C0504D', 'control': '#4F81BD', 'background_A': '#A8A8A8',
          'background_B': '#D8D8D8'}
ORDER = ['background_A', 'background_B', 'control', 'candidate']
LABEL = {'background_A': 'background, A copies', 'background_B': 'background, B copies',
         'control': 'family controls', 'candidate': 'candidates'}
SHORT = {'background_A': 'bg A', 'background_B': 'bg B', 'control': 'controls',
         'candidate': 'candidates'}
frac = {g: sum(1 for r in groups[g] if float(r['p']) < 0.05) / len(groups[g])
        for g in ORDER}

fig = plt.figure(figsize=(3.35, 3.1))          # half the former width, same height
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.62,
                      left=0.115, right=0.985, top=0.895, bottom=0.215)

# a - branch-model omega, with significance per set and the group rates inset
ax = fig.add_subplot(gs[0, 0])
lim = 1.5
for g in ORDER:
    bg_s = g.startswith('background')
    for r in groups[g]:
        x, y = float(r['omega_bg']), float(r['omega_rac'])
        if float(r['p']) < 0.05:                       # significant: filled
            ax.scatter([x], [y], s=9 if bg_s else 15, c=colors[g],
                       alpha=0.55 if bg_s else 0.9, lw=0, zorder=3)
        else:                                          # not significant: open
            ax.scatter([x], [y], s=9 if bg_s else 15, facecolor='white',
                       edgecolor=colors[g], lw=0.7, alpha=0.9, zorder=3)
    mx = np.mean([float(r['omega_bg']) for r in groups[g]])
    my = np.mean([float(r['omega_rac']) for r in groups[g]])
    ax.scatter([mx], [my], s=46, marker='D', c=colors[g], edgecolor='white',
               lw=0.8, zorder=5, label=f'{SHORT[g]} ({frac[g]:.0%})')
ax.plot([0, lim], [0, lim], ls='--', lw=0.7, color='#888888', zorder=2)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('background \u03c9', fontsize=6.6)
ax.set_ylabel('racemose-branch \u03c9', fontsize=6.6)
ax.tick_params(labelsize=5.8)
ax.set_xticks([0, 0.5, 1.0, 1.5])
ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.legend(fontsize=4.8, loc='upper left', handletextpad=0.3, labelspacing=0.30,
          borderpad=0.2)
ax.set_title('a  branch-model \u03c9', fontsize=7.4,
             loc='left')

axi = ax.inset_axes([0.50, 0.145, 0.48, 0.28])
axi.barh(np.arange(4), [frac[g] for g in ORDER],
         color=[colors[g] for g in ORDER], height=0.62, zorder=3)
axi.axvline(0.05, color='#C0504D', ls=':', lw=0.7, zorder=4)
for i, g in enumerate(ORDER):
    axi.text(frac[g] + 0.025, i, f'{frac[g]:.0%}', va='center', fontsize=4.2,
             color='#555555', zorder=5)
axi.set_xlim(0, 0.62)
axi.set_ylim(-0.6, 3.9)
axi.set_xticks([])
axi.set_yticks([])
axi.tick_params(labelsize=4.2, length=0)
axi.set_title('fraction with p < 0.05', fontsize=4.6, pad=1.2, color='#555555')
for sp in ('top', 'right'):
    axi.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    axi.spines[sp].set_linewidth(0.6)
axi.patch.set_alpha(0.92)

# b - ratio distribution
ax3 = fig.add_subplot(gs[0, 1])
for g, dy in [('background_A', -0.12), ('background_B', 0.12)]:
    vals = [np.log2(float(r['omega_rac']) / float(r['omega_bg']))
            for r in groups[g] if float(r['omega_bg']) > 0]
    ax3.scatter(vals, np.random.default_rng(1).normal(dy, 0.03, len(vals)),
                s=8, c=colors[g], alpha=0.7, lw=0)
cand_vals = [np.log2(float(r['omega_rac']) / float(r['omega_bg']))
             for r in groups['candidate'] if float(r['omega_bg']) > 0]
ax3.scatter(cand_vals, np.zeros(len(cand_vals)), s=15, c=colors['candidate'],
            alpha=0.9, lw=0)
ax3.axvline(0, color='#888888', ls='--', lw=0.7)
ax3.set_yticks([-0.12, 0.0, 0.12], ['bg A', 'cand.', 'bg B'], fontsize=5.0)
ax3.set_ylim(-0.42, 0.38)
ax3.set_xlabel('log2 \u03c9 ratio', fontsize=6.6)
ax3.tick_params(axis='x', labelsize=5.6)
ax3.set_title('b  per-set \u03c9 ratio', fontsize=7.4, loc='left')

fig.text(0.115, 0.018,
         'a: filled, p < 0.05; open, not significant;\ndiamonds, group means. '
         'b: one point per gene set.',
         fontsize=4.4, color='#555555', linespacing=1.5)

for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
    fig.savefig(FIG / f'fig6_dnds_branch.{ext}', bbox_inches='tight', **kw)
fig.savefig('${HOME}/AppData/Local/Temp/opencode/fig1e_preview/fig6_dnds_preview.png',
            dpi=150, bbox_inches='tight')
print('saved fig6 and', RES / 'branch_dnds_integrated.tsv')
