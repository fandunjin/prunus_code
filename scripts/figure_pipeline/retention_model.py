#!/usr/bin/env python3
"""retention_model.py - one schematic of the racemose-inflorescence retention model.

Core conclusion the figure defends:
  Inflorescence polarity matches genomic polarity in Prunus.
  raceme          = ancestral programme + polyploid buffering after one shared A x B WGD
  solitary/corymb = derived reduction on the contraction-biased B-associated diploid branch

Evidence encoded (all from the all-17-track analysis package):
  WGD stem branch   +11,747 / -95 families; 956 / 0 branch-significant (p < 0.05)
  464 families      retained above the WGD expectation (p < 0.05, >1.3x, >=2/3 species)
  retained modules  nucleus/chromatin, Golgi-vesicle, meiosis/chiasma, RNA processing
  regulators        BLH1, YABBY1/2/3, SVP/DAM, SPL9/15
  A/B asymmetry     A retains less than B; relaxed selection A > B
                    (44% vs 19% of random single-copy genes significant; 5% null)
  diploid branch    contraction-biased: 548/2,263 - 477/1,456 - 1,124/4,314 - 1,455/3,925

Output: retention_model.{svg,pdf,png}  (SVG/PDF keep live text for Illustrator)
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import (Arc, Circle, Ellipse, PathPatch, Polygon, Rectangle)

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

# ------------------------------------------------------------------ palette
C_A = '#C0504D'        # A subgenome
C_B = '#4F81BD'        # B subgenome
C_COR = '#5B8FD6'      # corymbose
C_SOL = '#7BAA5B'      # solitary
C_RAC = '#B03A2E'      # racemose
C_HYB = '#7B5EA7'      # WGD / hybrid
C_TXT = '#3C3C3C'
C_SUB = '#6E6E6E'
C_BOX = '#F4F4F2'

W, H = 183.0, 128.0
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')


# ------------------------------------------------------------------ helpers
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


def node(x, y, color=C_TXT, r=1.0, z=8):
    ax.add_patch(Circle((x, y), r, fc=color, ec='white', lw=0.6, zorder=z))


def star(x, y, r, color, rot=90, z=14, ec='white', lw=1.0, ratio=0.42):
    pts = []
    for i in range(10):
        rr = r if i % 2 == 0 else r * ratio
        a = math.radians(rot + i * 36)
        pts.append((x + rr * math.cos(a), y + rr * math.sin(a)))
    ax.add_patch(Polygon(pts, closed=True, fc=color, ec=ec, lw=lw, zorder=z))


def card(x, y, w, h, ec='#DDDDD8', fc=C_BOX, z=2, lw=0.7):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z))


def cable(x0, x1, y, z=6, lw_a=2.7, lw_b=1.3, off=0.9, nwaves=5.0):
    """A+B subgenome cable: a red A strand and a blue B strand braided together."""
    xs = [x0 + (x1 - x0) * i / 40.0 for i in range(41)]
    for sgn in (+1, -1):
        ys = [y + sgn * off * math.sin((x - x0) / (x1 - x0) * math.pi * nwaves)
              for x in xs]
        ax.plot(xs, ys, color=C_B if sgn > 0 else C_A, lw=lw_b,
                solid_capstyle='round', zorder=z)


def inflo_icon(x, y, kind, w=10.0, color=None, z=9):
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


def glyph_nucleus(x, y, s=1.0, col=C_A, z=10):
    ax.add_patch(Circle((x, y), 2.3 * s, fc='none', ec=col, lw=1.0, zorder=z))
    ax.add_patch(Circle((x + 0.6 * s, y + 0.3 * s), 0.85 * s, fc=col, ec='none',
                        alpha=0.55, zorder=z + 1))
    for i in range(6):
        a = math.radians(35 * i)
        ax.plot([x + 0.8 * s * math.cos(a), x + 1.9 * s * math.cos(a)],
                [y + 0.8 * s * math.sin(a), y + 1.9 * s * math.sin(a)],
                color=col, lw=0.5, alpha=0.7, zorder=z)


def glyph_golgi(x, y, s=1.0, col=C_B, z=10):
    for i, rr in enumerate((2.9, 2.1, 1.3)):
        ax.add_patch(Arc((x, y + i * 1.05 * s), rr * 2 * s, rr * 1.5 * s,
                         theta1=200, theta2=340, color=col, lw=1.05, zorder=z))
    for i in range(3):
        ax.add_patch(Circle((x - 2.0 * s + i * 2.0 * s, y - 2.1 * s), 0.5 * s,
                            fc=col, ec='none', alpha=0.8, zorder=z))


def glyph_chiasma(x, y, s=1.0, col=C_SOL, z=10):
    ax.plot([x - 2.4 * s, x + 2.4 * s], [y - 1.9 * s, y + 1.9 * s], color=col,
            lw=1.05, solid_capstyle='round', zorder=z)
    ax.plot([x - 2.4 * s, x + 2.4 * s], [y + 1.9 * s, y - 1.9 * s], color=col,
            lw=1.05, solid_capstyle='round', zorder=z)
    ax.add_patch(Circle((x, y), 0.5 * s, fc=col, ec='white', lw=0.4, zorder=z + 1))


# ================================================================== 1. history
node(13.0, 96.0, C_TXT, r=1.15)
inflo_icon(13.0, 104.0, 'raceme', w=9.5)
txt(13.0, 110.5, 'Prunus ancestor\nraceme  |  x = 9', fs=4.7, ha='center',
    va='bottom', color=C_TXT, lsp=1.25)

# A and B parental lineages (both converge on the shared WGD)
branch([(13.0, 96.0), (20.0, 102.5), (38.0, 102.5), (50.0, 99.0)], C_A, lw=2.6)
branch([(13.0, 96.0), (23.0, 84.0), (38.0, 84.0)], C_B, lw=2.6)
branch([(38.0, 84.0), (46.0, 90.0), (50.0, 94.5)], C_B, lw=2.6)
txt(25.0, 104.6, 'A lineage', fs=4.9, color=C_A, ha='center', weight='bold')
txt(25.0, 80.4, 'B lineage', fs=4.9, color=C_B, ha='center', weight='bold')

star(53.5, 96.5, 3.4, C_HYB)
txt(56.5, 102.8, 'shared WGD  (A \u00d7 B allopolyploidisation)\n> 42 Ma',
    fs=5.0, color=C_HYB, ha='center', va='bottom', weight='bold', lsp=1.35)

# racemose clade: A+B cable with three tips
cable(57.5, 96.0, 96.5, lw_a=2.7, lw_b=1.3, off=0.9)
for xt, nm in ((104.0, 'P. padus'), (124.0, 'P. spinulosa'), (144.0, 'P. zippeliana')):
    branch([(96.0, 96.5), (xt - 2.0, 96.5), (xt, 102.0)], C_RAC, lw=1.5)
    node(xt, 102.0, C_RAC, r=0.8)
    inflo_icon(xt, 110.0, 'raceme', w=9.0)
    txt(xt, 116.5, nm, fs=4.7, color=C_RAC, ha='center', va='bottom',
        style='italic')
txt(76.0, 91.5, 'racemose clade  (2n = 4x = 32)', fs=4.8, color=C_RAC,
    ha='center', va='top')

# B-associated diploid clade
branch([(38.0, 84.0), (56.0, 78.0), (76.0, 78.0)], C_B, lw=2.4)
branch([(76.0, 78.0), (86.0, 73.0)], C_B, lw=2.2)
branch([(76.0, 78.0), (108.0, 73.0)], C_B, lw=2.2)
for xt, kind, nm, col in ((86.0, 'corymb', 'corymbose', C_COR),
                          (108.0, 'solitary', 'solitary', C_SOL)):
    node(xt, 73.0, C_B, r=0.8)
    inflo_icon(xt, 65.0, kind, w=9.0)
    txt(xt, 59.0, nm, fs=4.9, ha='center', va='top', weight='bold', color=col)
txt(60.0, 76.0, 'diploid Prunus\n(2n = 2x = 16)', fs=4.6, color=C_B,
    ha='center', va='top', lsp=1.25)

# ============================================================ 2. genomic readout
YT, HT = 30.0, 24.0
card(4.0, YT, 50.0, HT)
txt(6.5, YT + HT - 2.6, 'WGD stem branch', fs=5.1, weight='bold', color=C_HYB)
txt(6.5, YT + HT - 8.0, '+11,747 / \u221295 families', fs=5.9, weight='bold',
    color=C_RAC)
txt(6.5, YT + HT - 12.2, '956 / 0 branch-significant (p < 0.05)', fs=4.4,
    color=C_SUB)
txt(6.5, YT + HT - 17.2, '464 families retained above the\nWGD expectation',
    fs=4.8, color=C_TXT, lsp=1.3)

card(57.0, YT, 66.0, HT)
txt(59.5, YT + HT - 2.6, 'Retained modules', fs=5.1, weight='bold', color=C_TXT)
rows = [(46.6, glyph_chiasma, C_SOL, 'meiosis  |  chiasma assembly',
         'ASY1  SPO11-1  PRD1  DMC1  RPA1  pif1'),
        (40.6, glyph_nucleus, C_A, 'nucleus  |  chromatin',
         'DNA-binding transcription repression'),
        (34.6, glyph_golgi, C_B, 'Golgi  |  vesicle transport',
         'ER-to-Golgi  |  trans-Golgi network')]
for gy, glyph, col, line1, line2 in rows:
    glyph(63.0, gy, s=0.95, col=col)
    txt(69.0, gy + 1.1, line1, fs=4.8, color=C_TXT)
    txt(69.0, gy - 1.9, line2, fs=4.1, color=C_SUB)

card(126.0, YT, 53.0, HT)
txt(128.5, YT + HT - 2.6, 'Meristem \u2013 floral regulators', fs=5.1,
    weight='bold', color=C_TXT)
txt(128.5, YT + HT - 8.2, 'BLH1  \u00b7  YABBY1/2/3', fs=5.5, weight='bold',
    color=C_RAC)
txt(128.5, YT + HT - 13.2, 'SVP/DAM  \u00b7  SPL9/15', fs=5.5, weight='bold',
    color=C_RAC)
ax.plot([128.5, 176.5], [YT + 8.0, YT + 8.0], color='#DDDDD8', lw=0.6, zorder=3)
txt(128.5, YT + 5.2, 'A subgenome retains less than B', fs=4.5, color=C_TXT)
txt(128.5, YT + 2.0, 'relaxed selection A > B  (44% vs 19%; null 5%)', fs=4.2,
    color=C_SUB)

# ================================================================ 3. polarity
card(4.0, 4.0, 175.0, 21.0, fc='white', ec='#DDDDD8')
inflo_icon(13.0, 14.5, 'raceme', w=10.5)
txt(23.5, 14.6, 'raceme', fs=7.2, weight='bold', color=C_RAC, va='center')
txt(40.0, 14.6, '=  ancestral inflorescence programme\n+  polyploid buffering',
    fs=5.4, color=C_TXT, va='center', lsp=1.35)

ax.plot([79.0, 79.0], [6.5, 22.5], color='#DDDDD8', lw=0.8, ls=(0, (2, 2)),
        zorder=3)

inflo_icon(86.0, 14.5, 'corymb', w=9.5)
inflo_icon(96.0, 14.5, 'solitary', w=9.5)
txt(105.0, 14.6, 'solitary / corymbose', fs=6.4, weight='bold', color=C_SOL,
    va='center')
txt(134.0, 16.6, '=  derived reduction on the contraction-\nbiased diploid branch',
    fs=5.4, color=C_TXT, va='center', lsp=1.35)
txt(134.0, 10.4, '548/2,263  \u00b7  477/1,456  \u00b7  1,124/4,314  \u00b7  1,455/3,925',
    fs=4.3, color=C_SUB, va='center')

# ------------------------------------------------------------------ export
for ext, kw in (('svg', {}), ('pdf', {}), ('png', dict(dpi=600))):
    fig.savefig(FIG / f'retention_model.{ext}', **kw)
print('saved retention_model.svg / .pdf / .png')
