#!/usr/bin/env python3
"""plot_fig6b_dnds_permutation_control.py - permuted-foreground control panel for Fig. 6b.

Left: mean foreground/background omega ratio, observed racemose (A, B) versus the
16 permuted foregrounds. Right: fraction of sets with p < 0.05, observed versus
permuted. Outputs fig6b_dnds_permutation_control.{pdf,png,svg}.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

B = Path('${PROJECT_ROOT}')
OUT = B / 'analysis/dnds_branch_20260922/results'
FIG = B / '0903/figures/inflorescence_evolution_review'

OBSERVED = {'A': {'mean_ratio': 1.72, 'frac_sig': 0.44}, 'B': {'mean_ratio': 1.46, 'frac_sig': 0.19}}
COLORS = {'A': '#C0504D', 'B': '#4F81BD'}


def read(fn):
    with (OUT / fn).open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 7,
        'axes.linewidth': 0.75, 'axes.spines.right': False, 'axes.spines.top': False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(3.6, 2.0), gridspec_kw={'wspace': 0.42})
    for ax, key, ylabel in ((axes[0], 'mean_ratio', 'Mean omega ratio'),
                            (axes[1], 'frac_sig', 'Fraction significant')):
        for i, g in enumerate(('A', 'B')):
            rows = read(f'branch_dnds_permutation_control_{g}.tsv')
            vals = [float(r[key]) for r in rows]
            x = i + np.random.RandomState(1).uniform(-0.12, 0.12, len(vals))
            ax.scatter(x, vals, s=7, color='#B0B0B0', edgecolor='none', zorder=2)
            ax.hlines(np.median(vals), i - 0.22, i + 0.22, color='#666666', lw=0.9, zorder=3)
            obs = OBSERVED[g][key]
            ax.scatter([i], [obs], marker='D', s=16, color=COLORS[g], zorder=4)
            ax.annotate(f'{obs:.2f}' if key == 'ratio' else f'{obs:.0%}',
                        (i, obs), textcoords='offset points', xytext=(0, 5),
                        ha='center', fontsize=6.2, color=COLORS[g])
        ax.set_xticks([0, 1], ['A', 'B'])
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel(ylabel)
        ax.tick_params(length=3, width=0.7)
        if key == 'frac':
            ax.set_ylim(bottom=0)
    axes[0].set_title('Permuted foregrounds (grey) versus observed', fontsize=6.6, loc='left')
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.86)
    for ext, kw in (('pdf', {}), ('svg', {}), ('png', {'dpi': 300})):
        fig.savefig(FIG / f'fig6b_dnds_permutation_control.{ext}', bbox_inches='tight', **kw)
    print('saved fig6b_dnds_permutation_control.{pdf,png,svg}')


if __name__ == '__main__':
    main()
