#!/usr/bin/env python3
"""plot_inflorescence_literature.py - schematic of the Prunus inflorescence
evolution as described in the literature.

Sources (text extracted from the PDFs):
  Su et al. 2023 (Plant Divers.): three clades by inflorescence type -
    racemose (Padus, Laurocerasus, Maddenia, Pygeum), solitary-flower
    (Amygdalus, Prunus, Microcerasus) and corymbose (Cerasus), the latter two
    sister; MRCA of Prunus racemose; solitary derived by reduction of flower
    number, corymb by suppression of the rachis; reversal to raceme in
    P. discadenia; crown ~67.3 Ma, solitary+corymbose split ~58.4 Ma.
  Wang et al. 2026 (Mol. Phylogenet. Evol.): crown ~63.0 Ma; racemose clade
    earliest-diverging (~61.4 Ma); corymbose + solitary split ~51.0 Ma; MRCA
    racemose and polyploid.
  Zhao et al. 2016 (PLoS ONE): three groups by inflorescence structure and
    ploidy - diploid solitary, diploid corymbose, polyploid racemose.

Output: inflorescence_evolution_literature.{svg,pdf,png}
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, PathPatch,
                                Polygon)

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'font.size': 6,
})

_CANDIDATES = [Path(r'${PROJECT_ROOT}'),
               Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'

C_RAC = '#B03A2E'
C_COR = '#5B8FD6'
C_SOL = '#7BAA5B'
C_TXT = '#3C3C3C'
C_SUB = '#6E6E6E'
C_STEM = '#555555'

W, H = 183.0, 122.0
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')


def txt(x, y, s, fs=5.0, color=C_TXT, ha='left', va='center', style='normal',
        weight='normal', z=30):
    return ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, style=style,
                   weight=weight, zorder=z)


def catmull(pts):
    P = [pts[0]] + list(pts) + [pts[-1]]
    verts, codes = [P[1]], [MPath.MOVETO]
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        c1 = (p1[0] + (p2[0] - p0[0]) / 6.0, p1[1] + (p2[1] - p0[1]) / 6.0)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6.0, p2[1] - (p3[1] - p1[1]) / 6.0)
        verts += [c1, c2, p2]
        codes += [MPath.CURVE4] * 3
    return MPath(verts, codes)


def branch(pts, color, lw=2.0, z=5, ls='-'):
    ax.add_patch(PathPatch(catmull(pts), fc='none', ec=color, lw=lw, ls=ls,
                           capstyle='round', joinstyle='round', zorder=z))


def curve_arrow(pts, color, lw=1.4, ms=9, z=7, ls='-'):
    ax.add_patch(FancyArrowPatch(path=catmull(pts), arrowstyle='-|>',
                                 mutation_scale=ms, color=color, lw=lw, zorder=z,
                                 ls=ls))


def node(x, y, color=C_TXT, r=1.0, z=8):
    ax.add_patch(Circle((x, y), r, fc=color, ec='white', lw=0.6, zorder=z))


def rod(x, y, h, color, lw=2.2, z=8):
    ax.plot([x, x], [y - h / 2, y + h / 2], color=color, lw=lw,
            solid_capstyle='round', zorder=z)


def chrom_row(x, y, n, color, gap=0.85, h=3.4, lw=2.0, z=8):
    for i in range(n):
        rod(x + i * gap, y, h, color, lw=lw, z=z)


def flower_icon(x, y, r, color='#F2A7BC', edge='#C97A93', center='#9CBF6B',
                z=9, dots='#B23A5A'):
    for i in range(5):
        a = math.radians(90 + i * 72)
        ax.add_patch(Circle((x + 0.58 * r * math.cos(a),
                             y + 0.58 * r * math.sin(a)), 0.46 * r, fc=color,
                            ec=edge, lw=0.4, zorder=z))
    ax.add_patch(Circle((x, y), 0.24 * r, fc=center, ec='none', zorder=z + 1))


def inflo_icon(x, y, kind, w=11.0, color=None, z=9):
    s = w / 10.0
    if kind == 'raceme':
        col = color or C_RAC
        ax.plot([x, x], [y - 4.4 * s, y + 3.8 * s], color='#6E8F4A', lw=2.4 * s,
                solid_capstyle='round', zorder=z)
        for i, yy in enumerate([-2.8, -1.7, -0.6, 0.5, 1.6, 2.7]):
            side = 1 if i % 2 == 0 else -1
            ax.add_patch(Circle((x + side * 1.7 * s, y + yy * s), 0.95 * s,
                                fc=col, ec='white', lw=0.4, zorder=z))
        ax.add_patch(Ellipse((x, y + 4.5 * s), 1.3 * s, 0.8 * s, fc='#8FBF6B',
                             ec='#5E8A45', lw=0.5, zorder=z))
    elif kind == 'corymb':
        col = color or C_COR
        ax.plot([x, x], [y - 4.4 * s, y - 0.8 * s], color='#8A6A4F', lw=2.6 * s,
                solid_capstyle='round', zorder=z)
        for dx in (-2.8, -1.0, 1.0, 2.8):
            ax.plot([x, x + dx * s], [y - 0.9 * s, y + 1.7 * s], color='#6E8F4A',
                    lw=1.5 * s, zorder=z)
            ax.add_patch(Circle((x + dx * s, y + 2.5 * s), 1.1 * s, fc=col,
                                ec='white', lw=0.4, zorder=z))
    elif kind == 'solitary':
        col = color or C_SOL
        ax.plot([x, x], [y - 4.4 * s, y - 0.6 * s], color='#6E8F4A', lw=2.4 * s,
                solid_capstyle='round', zorder=z)
        ax.add_patch(Circle((x, y + 1.7 * s), 2.3 * s, fc=col, ec='white',
                            lw=0.45, zorder=z))
        ax.add_patch(Circle((x, y + 1.7 * s), 0.9 * s, fc='#F2E3B3',
                            ec='#C9A227', lw=0.4, zorder=z + 1))
    elif kind == 'mini_raceme':          # 3 flowers on a short axis (for glyphs)
        col = color or C_RAC
        ax.plot([x, x], [y - 3.0 * s, y + 3.0 * s], color='#6E8F4A', lw=2.0 * s,
                solid_capstyle='round', zorder=z)
        for i, yy in enumerate([-1.9, 0.0, 1.9]):
            side = 1 if i % 2 == 0 else -1
            ax.add_patch(Circle((x + side * 1.5 * s, y + yy * s), 0.95 * s,
                                fc=col, ec='white', lw=0.4, zorder=z))
        flower_icon(x + 1.5 * s, y + 1.9 * s, 1.5 * s, z=z + 2)


# ================================================================== tree
# root: ancestral raceme
branch([(16.0, 58.0), (25.0, 58.0), (34.0, 58.0)], C_STEM, lw=2.6)
node(16.0, 58.0, color='#333333', r=1.3)
inflo_icon(16.0, 70.0, 'raceme', w=12.0, color=C_RAC)
txt(16.0, 46.5, '67.3 / 63.0 Ma', fs=4.6, color=C_SUB, ha='center')
node(34.0, 58.0, color='#333333', r=1.1)

# racemose clade (earliest-diverging, polyploid)
branch([(34.0, 58.0), (42.0, 66.0), (50.0, 80.0), (56.0, 88.0)], C_RAC, lw=2.6)
node(56.0, 88.0, color=C_RAC, r=1.1)
tips = [('Padus', 78.0), ('Laurocerasus', 86.0), ('Maddenia', 94.0),
        ('Pygeum', 102.0)]
for i, (name, yy) in enumerate(tips):
    branch([(56.0, 88.0), (64.0, 84.0 + i * 3.0), (76.0, yy), (90.0, yy)],
           C_RAC, lw=1.8)
    inflo_icon(94.5, yy + 0.2, 'raceme', w=7.5, color=C_RAC)
    txt(100.0, yy, name, fs=5.0, color=C_RAC, weight='bold')
chrom_row(36.0, 92.5, 8, C_RAC)
chrom_row(36.0, 89.0, 8, C_RAC)
txt(36.0, 85.0, '2n = 32', fs=4.6, color=C_RAC)

# (solitary + corymbose) stem
branch([(34.0, 58.0), (40.0, 55.0), (44.0, 52.0)], C_STEM, lw=2.6)
node(44.0, 52.0, color='#333333', r=1.1)
txt(44.0, 47.0, '58.4 / 51.0 Ma', fs=4.6, color=C_SUB, ha='center')

# corymbose clade (diploid)
branch([(44.0, 52.0), (50.0, 58.0), (58.0, 66.0)], C_COR, lw=2.4)
node(58.0, 66.0, color=C_COR, r=1.1)
for yy in (62.0, 70.0):
    branch([(58.0, 66.0), (70.0, yy), (90.0, yy)], C_COR, lw=1.8)
    inflo_icon(94.5, yy + 0.2, 'corymb', w=7.5, color=C_COR)
txt(100.0, 62.0, 'Cerasus', fs=5.0, color=C_COR, weight='bold')
txt(100.0, 70.0, 'P. maackii', fs=5.0, color=C_COR, weight='bold')
chrom_row(52.0, 74.0, 8, C_COR)
txt(52.0, 70.5, '2n = 16', fs=4.6, color=C_COR)

# solitary clade (diploid)
branch([(44.0, 52.0), (50.0, 42.0), (58.0, 34.0)], C_SOL, lw=2.4)
node(58.0, 34.0, color=C_SOL, r=1.1)
for yy, name in ((26.0, 'Amygdalus'), (34.0, 'Prunus'), (42.0, 'Microcerasus')):
    branch([(58.0, 34.0), (70.0, yy), (90.0, yy)], C_SOL, lw=1.8)
    inflo_icon(94.5, yy + 0.2, 'solitary', w=7.5, color=C_SOL)
    txt(100.0, yy, name, fs=5.0, color=C_SOL, weight='bold')
chrom_row(52.0, 24.0, 8, C_SOL)
txt(52.0, 20.5, '2n = 16', fs=4.6, color=C_SOL)

# ================================================================== transitions
# raceme -> simple flower (reduction of flower number)
inflo_icon(32.0, 40.0, 'mini_raceme', w=8.0, color=C_RAC)
curve_arrow([(37.0, 40.0), (41.0, 40.0)], '#4F6B33', lw=1.3, ms=9)
flower_icon(44.5, 40.0, 2.6, z=9)
curve_arrow([(44.0, 42.5), (45.5, 49.0)], '#999999', lw=0.7, ms=6, ls=(0, (2, 2)))

# simple flower -> corymb (suppression of the rachis)
flower_icon(64.0, 80.0, 2.6, z=9)
curve_arrow([(68.0, 80.0), (71.5, 80.0)], '#4F6B33', lw=1.3, ms=9)
inflo_icon(76.0, 80.0, 'corymb', w=7.5, color=C_COR)
curve_arrow([(70.0, 77.0), (66.5, 69.0)], '#999999', lw=0.7, ms=6, ls=(0, (2, 2)))

# reversal: P. discadenia (raceme again, rachis elongation)
curve_arrow([(106.0, 62.0), (112.0, 62.0)], C_RAC, lw=1.2, ms=8, ls=(0, (2, 2)))
inflo_icon(117.0, 62.2, 'raceme', w=7.5, color=C_RAC)
txt(121.5, 62.0, 'P. discadenia', fs=4.8, color=C_RAC)

# ------------------------------------------------------------------ export
for ext, kw in (('svg', {}), ('pdf', {}), ('png', dict(dpi=600))):
    fig.savefig(FIG / f'inflorescence_evolution_literature.{ext}', **kw)
print('saved inflorescence_evolution_literature.svg / .pdf / .png')
