#!/usr/bin/env python3
"""plot_evolution_model.py - one integrated schematic of the racemose Prunus
clade and the five lines of evidence for a single shared WGD.

No frames and (almost) no text: the tree carries the story, and every piece of
evidence is drawn directly next to the branch it supports.

  scenario   : Prunus ancestor (x = 9, raceme) -> A/B lineages -> B to the
               diploid corymbose/solitary clades; A x B -> one shared WGD
               (five-pointed star) -> racemose clade as an A+B subgenome cable.
  1 parallel : two identical A/B subgenome trees (P. padus first in both).
  2 Ks       : A/B Ks modes with the older A ages (41.8 / 19.5 vs 34.2 / 11.8).
  3 genes    : gene-family duplications peak at the racemose node
               (3,126 > 2,458 > 1,475).
  4 POInT    : pillar grids, multiple origins (x) vs one shared root (tick).
  5 dotplots : the same chromosome-scale signature in all three species.

Output: evolution_model.{svg,pdf,png}  (SVG keeps live text for Illustrator)
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, PathPatch,
                                Polygon, Rectangle)

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

C_A = '#C0504D'
C_B = '#4F81BD'
C_COR = '#5B8FD6'
C_SOL = '#7BAA5B'
C_RAC = '#B03A2E'
C_TXT = '#3C3C3C'
C_SUB = '#6E6E6E'
C_HYB = '#7B5EA7'

W, H = 183.0, 140.0
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')


def txt(x, y, s, fs=5.0, color=C_TXT, ha='left', va='center', style='normal',
        weight='normal', lsp=1.35, z=30):
    return ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, style=style,
                   weight=weight, linespacing=lsp, zorder=z)


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


def curve_arrow(pts, color, lw=1.6, ms=11, z=7):
    ax.add_patch(FancyArrowPatch(path=catmull(pts), arrowstyle='-|>',
                                 mutation_scale=ms, color=color, lw=lw, zorder=z))


def node(x, y, color=C_TXT, r=1.0, z=8):
    ax.add_patch(Circle((x, y), r, fc=color, ec='white', lw=0.6, zorder=z))


def star(x, y, r, color, rot=90, z=14, ec='white', lw=1.1, ratio=0.42):
    pts = []
    for i in range(10):
        rr = r if i % 2 == 0 else r * ratio
        a = math.radians(rot + i * 36)
        pts.append((x + rr * math.cos(a), y + rr * math.sin(a)))
    ax.add_patch(Polygon(pts, closed=True, fc=color, ec=ec, lw=lw, zorder=z))


def badge(x, y, letter, color, r=2.2, fs=5.6):
    ax.add_patch(Circle((x, y), r, fc=color, ec='white', lw=0.8, zorder=13))
    txt(x, y - 0.1, letter, fs=fs, color='white', ha='center', va='center',
        weight='bold', z=14)


def rod(x, y, h, color, seg=None, seg_color=None, lw=2.4, z=8):
    ax.plot([x, x], [y - h / 2, y + h / 2], color=color, lw=lw,
            solid_capstyle='round', zorder=z)
    if seg:
        y0 = y - h / 2 + seg[0] * h
        y1 = y - h / 2 + seg[1] * h
        ax.plot([x, x], [y0, y1], color=seg_color, lw=lw + 0.8,
                solid_capstyle='butt', zorder=z + 1)


def chrom_set(x, y, n, color, h=4.6, gap=1.05, segs=None, lw=2.4, z=8):
    for i in range(n):
        seg = segs.get(i) if segs else None
        rod(x + i * gap, y, h, color, seg=seg,
            seg_color='#F2C14E' if seg else None, lw=lw, z=z)


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
    else:
        col = color or C_SOL
        ax.plot([x, x], [y - 4.4 * s, y - 0.6 * s], color='#6E8F4A', lw=2.4 * s,
                solid_capstyle='round', zorder=z)
        ax.add_patch(Circle((x, y + 1.7 * s), 2.3 * s, fc=col, ec='white',
                            lw=0.45, zorder=z))
        ax.add_patch(Circle((x, y + 1.7 * s), 0.9 * s, fc='#F2E3B3',
                            ec='#C9A227', lw=0.4, zorder=z + 1))


def gauss(x0, x1, mu, sig, amp, color, z=6, base=0.0):
    xs = [x0 + (x1 - x0) * i / 60.0 for i in range(61)]
    ys = [base + amp * math.exp(-0.5 * ((x - mu) / sig) ** 2) for x in xs]
    ax.fill_between(xs, base, ys, color=color, alpha=0.28, lw=0, zorder=z)
    ax.plot(xs, ys, color=color, lw=1.1, zorder=z + 1)


def mini_tree(x0, y0, color, lw=1.4):
    branch([(x0, y0), (x0 + 6, y0)], color, lw=lw)
    branch([(x0 + 6, y0), (x0 + 11, y0 + 6), (x0 + 19, y0 + 6)], color, lw=lw)
    branch([(x0 + 6, y0), (x0 + 11, y0), (x0 + 11, y0 - 4)], color, lw=lw)
    branch([(x0 + 11, y0 - 4), (x0 + 15, y0 - 1), (x0 + 19, y0 - 1)], color, lw=lw)
    branch([(x0 + 11, y0 - 4), (x0 + 15, y0 - 6), (x0 + 19, y0 - 6)], color, lw=lw)


def dotplot_glyph(x, y, size):
    ax.add_patch(Rectangle((x, y), size, size, fc='white', ec='#C9C9C9', lw=0.6,
                           zorder=5))
    n = 8
    for i in range(n):
        px = x + size * (0.07 + 0.86 * i / (n - 1))
        py = y + size * (0.07 + 0.86 * i / (n - 1))
        ax.add_patch(Circle((px, py), 0.20, fc='#8A8A8A', ec='none', zorder=6))
    for (bx, by, w, h) in ((0.10, 0.50, 0.20, 0.20), (0.56, 0.10, 0.20, 0.18)):
        for i in range(4):
            for j in range(4):
                px = x + size * (bx + w * i / 3.0)
                py = y + size * (by + h * j / 3.0)
                ax.add_patch(Circle((px, py), 0.17, fc=C_RAC, ec='none',
                                    alpha=0.85, zorder=6))


def pillar_grid(x, y, w, h, ncol=6, nrow=4, z=6):
    for i in range(ncol):
        for j in range(nrow):
            ax.add_patch(Circle((x + w * i / (ncol - 1), y + h * j / (nrow - 1)),
                                0.30, fc='#CFCFCF', ec='none', zorder=z))


def leader(pts):
    branch(pts, '#BBBBBB', lw=0.7, z=4, ls=(0, (2.2, 2.2)))


# ================================================================== tree
# ancestor
chrom_set(9.0, 72.0, 9, '#8A8A8A', h=6.4, gap=1.15, lw=2.8)
txt(20.5, 72.0, 'x = 9', fs=5.2, color=C_SUB)
inflo_icon(13.5, 58.0, 'raceme', w=11.0, color=C_RAC)

# stem + translocation
branch([(19.5, 72.0), (32.0, 72.0), (44.0, 72.0)], '#555555', lw=2.2)
node(44.0, 72.0, r=0.9)
rod(28.4, 78.0, 5.4, '#8A8A8A', seg=(0.45, 0.75), seg_color=C_A, lw=2.2)
rod(31.6, 78.0, 5.4, '#8A8A8A', seg=(0.25, 0.55), seg_color=C_B, lw=2.2)
ax.plot([27.6, 32.4], [75.6, 80.4], color='#555555', lw=0.8, zorder=10)
ax.plot([27.6, 32.4], [80.4, 75.6], color='#555555', lw=0.8, zorder=10)

# A / B lineages
branch([(44.0, 72.0), (54.0, 75.0), (66.0, 84.0), (78.0, 96.0), (88.0, 106.0)],
       C_A, lw=2.2)
branch([(44.0, 72.0), (52.0, 68.0), (58.0, 58.0), (64.0, 48.0), (70.0, 42.0)],
       C_B, lw=2.2)
badge(50.0, 78.0, 'A', C_A, r=2.6, fs=7.0)
badge(48.0, 64.0, 'B', C_B, r=2.6, fs=7.0)

# karyotype changes (with leaders to their branches)
chrom_set(30.0, 88.0, 9, '#B9B9B9', h=3.6, gap=0.95, lw=1.8)
curve_arrow([(39.6, 88.0), (42.4, 88.0)], C_A, lw=1.2, ms=8)
chrom_set(45.0, 88.0, 8, C_A, h=3.6, gap=0.95,
          segs={1: (0.30, 0.70), 5: (0.20, 0.55)}, lw=2.0)
leader([(49.0, 90.0), (56.0, 78.0)])
chrom_set(30.0, 54.0, 9, '#B9B9B9', h=3.6, gap=0.95, lw=1.8)
curve_arrow([(39.6, 54.0), (42.4, 54.0)], C_B, lw=1.2, ms=8)
chrom_set(45.0, 54.0, 8, C_B, h=3.6, gap=0.95,
          segs={3: (0.25, 0.60)}, lw=2.0)
leader([(50.0, 56.0), (57.0, 60.0)])

# WGD star + doubling
star(88.0, 106.0, 6.2, C_RAC)
chrom_set(70.0, 121.0, 8, C_A, h=3.2, gap=0.80, lw=1.8)
chrom_set(70.0, 117.5, 8, C_B, h=3.2, gap=0.80, lw=1.8)
curve_arrow([(77.2, 119.2), (80.2, 119.2)], C_RAC, lw=1.2, ms=8)
for i in range(16):
    rod(83.0 + i * 0.80, 119.2, 3.2, C_A if i % 2 == 0 else C_B, lw=1.8)
txt(97.0, 119.2, '2n = 4x = 32', fs=4.8, color=C_RAC)

# hybridisation arrow from the B parent
curve_arrow([(54.0, 62.0), (60.0, 76.0), (64.0, 90.0), (80.0, 101.0)], C_HYB,
            lw=1.8, ms=12)


# racemose clade as an A+B cable
def cable(pts, z=6, lw_a=2.9, lw_b=1.3, off=1.35):
    dn = [(p[0], p[1] - off) for p in pts]
    branch(pts, C_A, lw=lw_a, z=z)
    branch(dn, C_B, lw=lw_b, z=z + 1)


trunk = [(88.0, 106.0), (96.0, 106.0), (108.0, 106.0)]
padus = [(94.0, 106.0), (106.0, 110.0), (116.0, 116.0), (128.0, 116.0),
         (140.0, 116.0)]
spin = [(108.0, 106.0), (118.0, 110.0), (128.0, 111.0), (140.0, 111.0)]
zipp = [(108.0, 106.0), (118.0, 103.0), (128.0, 102.0), (140.0, 102.0)]
for pts in (trunk, padus, spin, zipp):
    cable(pts)
txt(92.0, 108.4, 'A', fs=5.0, color=C_A, weight='bold')
txt(92.0, 102.4, 'B', fs=5.0, color=C_B, weight='bold')
for yy, name in ((116.0, 'P. padus'), (111.0, 'P. spinulosa'),
                 (102.0, 'P. zippeliana')):
    inflo_icon(142.5, yy + 0.2, 'raceme', w=5.6, color=C_RAC)
    txt(146.5, yy, name, fs=5.4, color=C_RAC, weight='bold')

# diploid descendants of the B lineage
branch([(70.0, 42.0), (82.0, 46.0), (96.0, 48.0), (118.0, 50.0), (140.0, 50.0)],
       C_COR, lw=2.2)
branch([(70.0, 42.0), (82.0, 39.0), (96.0, 36.0), (118.0, 32.0), (140.0, 32.0)],
       C_SOL, lw=2.2)
node(70.0, 42.0, color='#555555', r=0.9)
inflo_icon(144.0, 50.2, 'corymb', w=7.0, color=C_COR)
inflo_icon(144.0, 32.2, 'solitary', w=7.0, color=C_SOL)
txt(131.0, 54.0, '2n = 16', fs=4.6, color=C_COR)
txt(131.0, 27.6, '2n = 16', fs=4.6, color=C_SOL)
txt(150.0, 52.0, 'Cerasus', fs=4.8, color=C_COR)
txt(150.0, 30.0, 'Armeniaca / Amygdalus', fs=4.8, color=C_SOL)
curve_arrow([(74.0, 44.0), (80.0, 47.0)], '#7BAA5B', lw=1.1, ms=8)
curve_arrow([(74.0, 40.0), (80.0, 37.0)], '#7BAA5B', lw=1.1, ms=8)

# ================================================================== evidence
# 1 parallel subgenome trees (right of the racemose tips)
mini_tree(150.0, 84.0, C_A)
mini_tree(150.0, 68.0, C_B)
ax.plot([161.0, 161.0], [80.0, 64.0], color='#C4C4C4', lw=0.7, ls=(0, (2, 2)),
        zorder=4)
ax.plot([156.0, 156.0], [84.0, 68.0], color='#C4C4C4', lw=0.7, ls=(0, (2, 2)),
        zorder=4)
txt(172.5, 90.0, '41.8', fs=4.6, color=C_A)
txt(172.5, 74.0, '34.2', fs=4.6, color=C_B)
txt(170.5, 90.0, 'P', fs=4.0, color=C_SUB)
txt(170.5, 83.0, 'S', fs=4.0, color=C_SUB)
txt(170.5, 78.0, 'Z', fs=4.0, color=C_SUB)
txt(170.5, 74.0 - 0.0, '', fs=4.0)
badge(147.5, 92.0, '1', '#9A9A9A')

# 2 Ks modes + ages (under the racemose cable)
gauss(106.0, 134.0, 116.0, 4.6, 9.0, C_B, z=6, base=86.0)
gauss(106.0, 134.0, 126.0, 4.6, 9.5, C_A, z=6, base=86.0)
txt(106.0, 82.5, 'Ks', fs=4.4, color=C_SUB)
for (yy, col, val, xmax) in [(78.0, C_A, '41.8', 21.0), (74.0, C_B, '34.2', 17.0),
                             (70.0, C_A, '19.5', 10.0), (66.0, C_B, '11.8', 6.0)]:
    ax.add_patch(Rectangle((106.0, yy - 0.8), xmax, 1.6, fc=col, ec='none',
                           zorder=6))
    txt(106.0 + xmax + 1.0, yy, val, fs=4.6, color=col, va='center')
leader([(130.0, 96.0), (128.0, 100.5)])
badge(104.0, 96.0, '2', '#9A9A9A')

# 3 gene-family duplication mapping (above the racemose clade)
branch([(100.0, 127.0), (106.0, 127.0)], '#777777', lw=1.4)
branch([(106.0, 127.0), (111.0, 131.0), (116.0, 134.0)], C_RAC, lw=1.8)
branch([(116.0, 134.0), (122.0, 138.0), (127.0, 138.0)], C_RAC, lw=1.8)
branch([(116.0, 134.0), (120.0, 134.0), (120.0, 130.0)], C_RAC, lw=1.8)
branch([(120.0, 130.0), (123.0, 132.0), (127.0, 132.0)], C_RAC, lw=1.8)
branch([(120.0, 130.0), (123.0, 127.0), (127.0, 127.0)], C_RAC, lw=1.8)
branch([(106.0, 127.0), (113.0, 122.0), (127.0, 122.0)], '#999999', lw=1.4)
for (cx, cy, r, lab) in ((106.0, 127.0, 3.1, '2,458'),
                         (116.0, 134.0, 3.6, '3,126'),
                         (120.0, 130.0, 2.3, '1,475')):
    ax.add_patch(Circle((cx, cy), r, fc='#FBE3D6', ec=C_RAC, lw=1.0, zorder=9))
    txt(cx, cy - 0.1, lab, fs=3.8, color='#7B2C29', ha='center', va='center',
        z=10)
txt(128.0, 138.0, 'P', fs=4.0, color=C_SUB)
txt(128.0, 132.0, 'S', fs=4.0, color=C_SUB)
txt(128.0, 127.0, 'Z', fs=4.0, color=C_SUB)
leader([(113.0, 124.0), (112.0, 118.0)])
badge(98.0, 129.0, '3', '#9A9A9A')

# 4 POInT model comparison (left of the star)
pillar_grid(32.0, 126.0, 15.0, 9.0)
for (sx, sy) in ((34.0, 133.0), (39.0, 129.5), (45.0, 134.0)):
    star(sx, sy, 1.7, '#B9B9B9', ec='white', lw=0.6, z=10)
ax.text(50.0, 130.5, '\u2717', fontsize=9.0, color='#C0392B', ha='center',
        va='center', zorder=20)
pillar_grid(32.0, 112.0, 15.0, 9.0)
star(39.5, 116.5, 2.4, C_RAC, ec='white', lw=0.8, z=12)
ax.text(50.0, 116.5, '\u2713', fontsize=9.0, color='#2E8B57', ha='center',
        va='center', zorder=20)
txt(56.0, 112.0, '\u0394AIC 565', fs=4.4, color=C_RAC)
txt(30.0, 137.5, '13,225', fs=4.4, color=C_SUB)
leader([(62.0, 124.0), (78.0, 112.0)])
badge(28.0, 135.0, '4', '#9A9A9A')

# 5 dot plots (same signature at every racemose tip)
dotplot_glyph(162.0, 112.5, 7.0)
dotplot_glyph(162.0, 107.5, 7.0)
dotplot_glyph(162.0, 98.5, 7.0)
ax.plot([160.5, 160.5], [98.5, 119.5], color='#C4C4C4', lw=0.7, ls=(0, (2, 2)),
        zorder=4)
badge(159.0, 122.0, '5', '#9A9A9A')

# ------------------------------------------------------------------ export
for ext, kw in (('svg', {}), ('pdf', {}), ('png', dict(dpi=600))):
    fig.savefig(FIG / f'evolution_model.{ext}', **kw)
print('saved evolution_model.svg / .pdf / .png')
