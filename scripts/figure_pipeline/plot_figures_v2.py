#!/usr/bin/env python3
"""plot_figures_v2.py - annotated figure set for the racemose-specific scan.

fig1  CAFE5 branch landscape (rebuilt from the 0903-based OrthoFinder run)
fig2  genome-wide scan volcano with functional categories
fig3  top retained families with a category strip
fig4  GO enrichment
fig5  corrected synthesis model (WGD on the racemose stem only; PCAM/PPE not affected)
fig0  integrated summary figure (a-f)
"""
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle
from scipy.stats import binomtest

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 7,
    'axes.spines.right': False, 'axes.spines.top': False, 'axes.linewidth': 0.8,
    'legend.frameon': False,
})

B = Path('${PROJECT_ROOT}')
SCAN = B / 'analysis/racemose_inflorescence_scan_20260922'
CAFE_NEW = B / 'analysis/cafe5/main_sep14/k3/Gamma_clade_results.txt'
CAFE_OLD = B / '0903/manuscript_pipeline/stage2_write/phase4_draft/figure5/output/source_data/cafe5_clade_branch_summary.tsv'
FIG = B / '0903/figures/inflorescence_evolution_review'
PREV = Path('${HOME}/AppData/Local/Temp/opencode/fig1e_preview')
INFL = PREV  # cropped inflorescence drawings

TOT = {'RAC': 186830, 'COR': 131314, 'SOL': 150034, 'OUT': 56748}

CATS = [
    ('TE / transposon', '#E8A33D',
     r'transposon|retrotransposon|retrovirus|transposase|integrase|reverse transcriptase|'
     r'pol polyprotein|gag|gypsy|copia|rnase h|retro'),
    ('floral / meristem / hormone', '#C0504D',
     r'flower|floral|infloresc|meristem|mads|knox|bel1|bell|svp|short vegetative|dam|'
     r'spl9|yabby|ap1|ap2|leafy|wus|clavata|tfl|ft-like|sepallata|soc1|agamous|pistillata|'
     r'auxin|cytokinin|gibberellin|abscisic|ethylene|jasmon|brassinosteroid|strigolactone|'
     r'photoperiod|circadian|dormancy|vernal|flowering'),
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


def classify(sym, desc):
    text = f'{sym} {desc}'.lower()
    for name, color, pat in CATS:
        if re.search(pat, text):
            return name, color
    return ('annotated (other)', '#BBBBBB') if (sym or desc) else ('unannotated', '#E6E6E6')


# ---------------------------------------------------------------- tree utilities
CAFE_TREE = B / 'analysis/cafe5/main_sep14/mcmctree_dated_ultrametric_tree_Ma.tre'
TAXON_NAME = {'pzi': 'P. zippeliana', 'psp': 'P. spinulosa', 'ppa': 'P. padus',
              'pcam': 'P. campanulata', 'ppe': 'P. persica',
              'gtr': 'G. trifoliata', 'rch': 'R. chinensis'}
TAXON_COLOR = {'pzi': '#C0504D', 'psp': '#C0504D', 'ppa': '#C0504D',
               'pcam': '#4F81BD', 'ppe': '#4F81BD', 'gtr': '#888888', 'rch': '#888888'}
INFL_TYPE = {'pzi': 'raceme', 'psp': 'raceme', 'ppa': 'raceme',
             'pcam': 'corymb', 'ppe': 'solitary'}


def parse_newick(s):
    s = s.strip().rstrip(';')
    pos = [0]

    def parse():
        if s[pos[0]] == '(':
            pos[0] += 1
            children = []
            while True:
                children.append(parse())
                if s[pos[0]] == ',':
                    pos[0] += 1
                else:
                    pos[0] += 1
                    break
            node = {'children': children}
        else:
            node = {'children': []}
        m = re.match(r'\s*([^(),:;]*)\s*(?::\s*([0-9.eE+-]+))?', s[pos[0]:])
        node['name'] = m.group(1).strip()
        node['len'] = float(m.group(2)) if m.group(2) else 0.0
        pos[0] += m.end()
        return node

    return parse()


def tree_prep(root):
    if not root['children']:
        root['h'] = 0.0
    else:
        for c in root['children']:
            tree_prep(c)
        root['h'] = max(c['h'] + c['len'] for c in root['children'])
    return root['h']


def tree_y(root, state):
    if not root['children']:
        root['y'] = state[0]
        state[0] += 1
    else:
        for c in root['children']:
            tree_y(c, state)
        root['y'] = sum(c['y'] for c in root['children']) / len(root['children'])


def draw_chronogram(ax, root, H, *, annotate=None, tags=None, xpad=22.0, lw=0.9,
                    fs_label=6.0, fs_annot=5.0, branch_style=None, show_axis=True):
    """Draw the dated tree with time increasing to the right (root at left, present at right).
    annotate: dict taxon_or_node -> (up, down, up_sig, down_sig)."""
    x_of = lambda h: H - h
    Y = lambda n: n['y']

    def color_lw(node):
        if branch_style:
            return branch_style(node)
        return '#555555', lw

    def rec(node):
        if node['children']:
            x = x_of(node['h'])
            ys = [Y(c) for c in node['children']]
            ax.plot([x, x], [min(ys), max(ys)], color='#555555', lw=0.8, zorder=1)
            for c in node['children']:
                col, w = color_lw(c)
                ax.plot([x, x_of(c['h'])], [Y(c), Y(c)], color=col, lw=w,
                        solid_capstyle='round', zorder=2)
                rec(c)

    rec(root)
    for node in iter_nodes(root):
        if not node['children']:
            nm = node['name']
            if nm:
                ax.text(x_of(0) + 1.5, Y(node), TAXON_NAME.get(nm, nm), fontsize=fs_label,
                        va='center', color=TAXON_COLOR.get(nm, '#333333'))
    if annotate:
        for node in iter_nodes(root):
            nm = node['name']
            key = nm if nm else f"<{node['id']}>"
            if key in annotate:
                u, d, us, ds = annotate[key]
                x0 = x_of(node['h'] + node['len'])
                x1 = x_of(node['h'])
                if not node['children']:
                    xm = x1 - 0.45 * (x1 - x0)
                    ax.text(xm, Y(node) + 0.58, f'+{u:,} ({us:,})', ha='center', va='bottom',
                            fontsize=fs_annot, color='#7B2C29')
                    ax.text(xm, Y(node) + 0.30, f'\u2212{d:,} ({ds:,})', ha='center', va='bottom',
                            fontsize=fs_annot, color='#2C4E7B')
                else:
                    xm = (x0 + x1) / 2
                    ax.text(xm, Y(node) + 0.28, f'+{u:,} ({us:,})', ha='center', va='bottom',
                            fontsize=fs_annot, color='#7B2C29')
                    ax.text(xm, Y(node) - 0.28, f'\u2212{d:,} ({ds:,})', ha='center', va='top',
                            fontsize=fs_annot, color='#2C4E7B')
    if tags:
        for node in iter_nodes(root):
            nm = node['name']
            key = nm if nm else f"<{node['id']}>"
            if key in tags:
                xm = ((H - (node['h'] + node['len'])) + (H - node['h'])) / 2
                txt, col = tags[key]
                ax.text(xm, Y(node), txt, ha='center', va='center', fontsize=5.0,
                        color=col, zorder=6,
                        bbox=dict(fc='white', ec='none', pad=0.6))
    if show_axis:
        ax.set_xlim(-2, H + xpad)
        ax.set_ylim(-0.8, 6.8)
        ax.set_yticks([])
        ticks = [100, 80, 60, 40, 20, 0]
        ax.set_xticks([x_of(t) for t in ticks], [str(t) for t in ticks])
        ax.set_xlabel('Time (Ma)', fontsize=7.0)
        ax.tick_params(axis='x', labelsize=6.4)
        for s in ('left', 'right', 'top'):
            ax.spines[s].set_visible(False)
    else:
        ax.set_xlim(-2, H + xpad)
        ax.set_ylim(-0.8, 6.8)
        ax.axis('off')


def iter_nodes(root):
    yield root
    for c in root['children']:
        yield from iter_nodes(c)


def id_nodes(root):
    """Assign CAFE5 node ids: tips 1..n in tree order, internal nodes n+1.. in preorder."""
    n_tips = [0]
    pre = [0]

    def visit(node):
        if not node['children']:
            n_tips[0] += 1
            node['id'] = n_tips[0]
        else:
            pre[0] += 1
            node['id'] = None
            for c in node['children']:
                visit(c)

    visit(root)
    counter = [n_tips[0]]
    node_ids = {}

    def assign(node):
        if node['children']:
            counter[0] += 1
            node['id'] = counter[0]
        node_ids[node['id']] = node
        for c in node['children']:
            assign(c)

    assign(root)
    return node_ids


def load_tree():
    root = parse_newick(CAFE_TREE.read_text())
    H = tree_prep(root)
    tree_y(root, [0])
    maxy = max(n['y'] for n in iter_nodes(root) if not n['children'])
    for n in iter_nodes(root):
        n['y'] = maxy - n['y']
    ids = id_nodes(root)
    return root, H, ids


# ---------------------------------------------------------------- inflorescence icons
def draw_flower(ax, cx, cy, r, petal='#F2A7BC', edge='#C97A93', center='#9CBF6B',
                dots='#B23A5A', n_pet=5, rot=0.0, lw=0.8, ndots=8):
    for i in range(n_pet):
        ang = rot + i * 2 * np.pi / n_pet
        px = cx + 0.55 * r * np.cos(ang)
        py = cy + 0.55 * r * np.sin(ang)
        ax.add_patch(Ellipse((px, py), 1.15 * r, 0.95 * r, angle=np.degrees(ang),
                             fc=petal, ec=edge, lw=lw, zorder=3))
    ax.add_patch(Circle((cx, cy), 0.22 * r, fc=center, ec='#6E8F4A', lw=0.5, zorder=4))
    for i in range(ndots):
        ang = i * 2 * np.pi / ndots
        ax.add_patch(Circle((cx + 0.55 * r * np.cos(ang), cy + 0.55 * r * np.sin(ang)),
                            0.07 * r, fc=dots, ec='none', zorder=4))


def draw_inflo(ax, kind, cx=0.5, cy=0.5, s=1.0):
    """Schematic inflorescence icon centred at (cx, cy) in axes coordinates."""
    stem = '#6E8F4A'; stem_dk = '#4F6B33'; wood = '#8A6A4F'
    if kind == 'raceme':
        x0 = cx
        ax.plot([x0, x0], [cy - 0.42 * s, cy + 0.40 * s], color=stem, lw=2.6 * s,
                solid_capstyle='round', zorder=2)
        ax.add_patch(Ellipse((x0, cy + 0.43 * s), 0.06 * s, 0.10 * s, fc='#C8B08A',
                             ec=wood, lw=0.6, zorder=3))
        ys = np.linspace(cy + 0.30 * s, cy - 0.30 * s, 6)
        for i, y in enumerate(ys):
            side = 1 if i % 2 == 0 else -1
            px = x0 + side * 0.11 * s
            ax.plot([x0, px], [y, y + 0.02 * s], color=stem, lw=1.4 * s, zorder=2)
            draw_flower(ax, px + side * 0.05 * s, y + 0.03 * s, 0.085 * s,
                        petal='#FDF6F8', edge='#C9BFA8', dots='#C9A227',
                        center='#9CBF6B', lw=0.7)
        for side in (-1, 1):
            ax.add_patch(Ellipse((x0 + side * 0.16 * s, cy - 0.40 * s), 0.30 * s, 0.10 * s,
                                 angle=side * -18, fc='#9CBF6B', ec=stem_dk, lw=0.7, zorder=2))
    elif kind == 'corymb':
        x0 = cx
        ax.plot([x0, x0], [cy - 0.42 * s, cy - 0.06 * s], color=wood, lw=2.8 * s,
                solid_capstyle='round', zorder=2)
        base = (x0, cy - 0.04 * s)
        tips = [(x0 - 0.30 * s, cy + 0.16 * s), (x0 - 0.11 * s, cy + 0.22 * s),
                (x0 + 0.10 * s, cy + 0.22 * s), (x0 + 0.30 * s, cy + 0.16 * s)]
        for tx, ty in tips:
            ax.plot([base[0], tx], [base[1], ty], color=stem, lw=1.3 * s, zorder=2)
        for tx, ty in tips:
            draw_flower(ax, tx, ty + 0.07 * s, 0.11 * s, petal='#F2A7BC',
                        edge='#C97A93', dots='#C9A227', center='#9CBF6B', lw=0.7)
        ax.add_patch(Ellipse((x0 - 0.10 * s, cy - 0.05 * s), 0.12 * s, 0.07 * s,
                             angle=25, fc='#9CBF6B', ec=stem_dk, lw=0.6, zorder=2))
        ax.add_patch(Ellipse((x0 + 0.09 * s, cy - 0.045 * s), 0.10 * s, 0.06 * s,
                             angle=-20, fc='#C8B08A', ec=wood, lw=0.6, zorder=2))
    elif kind == 'solitary':
        x0 = cx - 0.05 * s
        ax.plot([x0, x0], [cy - 0.40 * s, cy - 0.02 * s], color=stem, lw=2.6 * s,
                solid_capstyle='round', zorder=2)
        draw_flower(ax, x0, cy + 0.16 * s, 0.16 * s, petal='#F2A7BC', edge='#C97A93',
                    dots='#B23A5A', center='#9CBF6B', lw=0.9)
        ax.add_patch(Ellipse((x0 - 0.16 * s, cy - 0.16 * s), 0.075 * s, 0.12 * s,
                             angle=-30, fc='#F2A7BC', ec='#C97A93', lw=0.7, zorder=2))
        ax.add_patch(Ellipse((x0 - 0.05 * s, cy - 0.13 * s), 0.11 * s, 0.05 * s,
                             angle=15, fc='#C8B08A', ec=wood, lw=0.6, zorder=2))
        ax.add_patch(Ellipse((x0 + 0.09 * s, cy - 0.13 * s), 0.07 * s, 0.11 * s,
                             angle=30, fc='#9CBF6B', ec=stem_dk, lw=0.6, zorder=2))


def icon_at(ax, x, y, w_in, kind, dy_in=0.0):
    """Place a vector inflorescence icon centred at data coordinates (x, y); w_in in inches."""
    fig = ax.figure
    pos = ax.get_position()
    ax_w_in = pos.width * fig.get_figwidth()
    ax_h_in = pos.height * fig.get_figheight()
    w_ax = w_in / ax_w_in
    h_ax = w_in * 1.3 / ax_h_in
    xa, ya = ax.transData.transform((x, y))
    xf, yf = ax.transAxes.inverted().transform((xa, ya))
    ia = ax.inset_axes([xf - w_ax / 2, yf - h_ax / 2 + dy_in / ax_h_in, w_ax, h_ax])
    ia.set_xlim(0, 1)
    ia.set_ylim(0, 1)
    ia.axis('off')
    draw_inflo(ia, kind, 0.5, 0.5, 1.0)
    return ia


# ---------------------------------------------------------------- fig1
def load_cafe5():
    """Return dict branch -> (up, down), preferring the rebuilt Sep14 run."""
    if CAFE_NEW.exists():
        d = {}
        with CAFE_NEW.open() as fh:
            rd = csv.reader(fh, delimiter='\t')
            next(rd)
            for row in rd:
                if len(row) < 3:
                    continue
                d[row[0]] = (int(row[1]), int(row[2]))
        return d, 'Sep14 (0903)'
    d = {}
    with CAFE_OLD.open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            d[r['branch']] = (int(r['families_changed_up']), int(r['families_changed_down']))
    return d, 'Aug25 (old)'


def load_cafe5_sig():
    """Branch-specific significant changes: stars in Gamma_asr.tre (branch p<0.05),
    split into up/down by the sign of Gamma_change.tab."""
    kdir = CAFE_NEW.parent
    asr, changef = kdir / 'Gamma_asr.tre', kdir / 'Gamma_change.tab'
    if not (asr.exists() and changef.exists()):
        return {}
    id2key = {}
    with CAFE_NEW.open() as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 3 or row[0].startswith('#'):
                continue
            m = re.match(r'([a-z]+)<(\d+)>', row[0]) or re.match(r'<(\d+)>', row[0])
            id2key[m.group(m.lastindex)] = row[0]
    star_fams = defaultdict(set)
    token = re.compile(r'(?:[a-z]+)?<(\d+)>(\*?)_')
    with asr.open(errors='replace') as fh:
        for line in fh:
            if not line.strip().startswith('TREE'):
                continue
            fam = line.split('=')[0].split()[-1]
            for nid, st in token.findall(line):
                if st and nid in id2key:
                    star_fams[id2key[nid]].add(fam)
    change = {}
    with changef.open() as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            change[p[0]] = dict(zip(hdr[1:], map(int, p[1:])))
    out = {}
    for b, fams in star_fams.items():
        out[b] = (sum(1 for f in fams if change.get(f, {}).get(b, 0) > 0),
                  sum(1 for f in fams if change.get(f, {}).get(b, 0) < 0))
    return out


def _branch_keys(data):
    key_by_id = {}
    for k in data:
        m = re.match(r'([a-z]+)<(\d+)>', k)
        if m:
            key_by_id[m.group(2)] = k
        else:
            m2 = re.match(r'<(\d+)>', k)
            if m2:
                key_by_id[m2.group(1)] = k
    return key_by_id


def fig1(ax=None, compact=False):
    """Dated tree with per-branch family expansion/contraction counts (CAFE5)."""
    data, source = load_cafe5()
    sig = load_cafe5_sig()
    root, H, ids = load_tree()
    key_by_id = _branch_keys(data)

    annotate = {}
    for nid, node in ids.items():
        k = key_by_id.get(str(nid))
        if k and k in data:
            u, d = data[k]
            su, sd = sig.get(k, (0, 0))
            annotate[node['name'] if node['name'] else f'<{nid}>'] = (u, d, su, sd)

    def branch_style(node):
        nm = node['name']
        k = nm if nm else f"<{node['id']}>"
        if k in data:
            u, d = data[k]
            frac = u / max(u + d, 1)
            if k == '<11>':
                return '#B03A2E', 2.4
            if frac >= 0.6:
                return '#D98B80', 1.8
            if frac <= 0.4:
                return '#7FA6D4', 1.8
        return '#666666', 0.9

    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(5.6, 2.5))
    draw_chronogram(ax, root, H, annotate=annotate, tags={'<11>': ('WGD', '#B03A2E')},
                    branch_style=branch_style, xpad=24,
                    fs_label=6.2 if not compact else 5.4,
                    fs_annot=5.0 if not compact else 4.3)
    if own:
        handles = [Line2D([], [], color=c, lw=2, label=l) for c, l in
                   [('#B03A2E', 'racemose stem (WGD)'), ('#D98B80', 'expansion-biased'),
                    ('#7FA6D4', 'contraction-biased')]]
        ax.legend(handles=handles, fontsize=5.4, loc='upper center',
                  bbox_to_anchor=(0.5, -0.28), ncol=3, handletextpad=0.4)
        ax.set_title('Family-size change along the dated tree (CAFE5; +up (\u2212down), '
                     'parentheses = significant p<0.05)', fontsize=7.6, loc='left')
        save(fig, 'fig1_cafe5_branch_landscape')
    return ax


def fig1_bars(ax=None):
    """Supplementary bar-chart version of the CAFE5 branch landscape."""
    data, source = load_cafe5()
    sig = load_cafe5_sig()
    labelmap = {'<11>': 'racemose stem\n(WGD)', '<12>': 'psp–pzi stem',
                '<10>': 'stem (rac.+pcam/ppe)', '<13>': 'pcam–ppe stem',
                '<9>': 'gtr–Prunus stem'}
    order = ['<11>', '<12>', '<10>', '<13>', '<9>', 'pzi', 'psp', 'ppa',
             'pcam', 'ppe', 'gtr', 'rch']
    key = {}
    for k in data:
        m = re.match(r'([a-z]+)<(\d+)>', k)
        if m:
            key[m.group(1)] = k
        elif k.startswith('<'):
            key[k] = k
    labels, ups, downs, sups, sdowns, tcolors = [], [], [], [], [], []
    for b in order:
        k = key.get(b, b)
        if k not in data:
            continue
        u, dn = data[k]
        su, sd = sig.get(k, (0, 0))
        labels.append(labelmap.get(b, b.upper()))
        ups.append(u)
        downs.append(dn)
        sups.append(su)
        sdowns.append(sd)
        tcolors.append('#C0504D' if b in ('<11>', 'pzi', 'psp', 'ppa') else '#333333')
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(3.9, 3.0))
    y = np.arange(len(labels))[::-1]
    ax.barh(y, ups, color='#F1C9C4', height=0.62)
    ax.barh(y, [-d for d in downs], color='#CBDDF0', height=0.62)
    ax.barh(y, sups, color='#B03A2E', height=0.62)
    ax.barh(y, [-d for d in sdowns], color='#3B6EA5', height=0.62)
    ax.axvline(0, color='black', lw=0.6)
    ax.set_yticks(y, labels, fontsize=6.2)
    for tick, c in zip(ax.get_yticklabels(), tcolors):
        tick.set_color(c)
    ax.set_xlabel('CAFE5 families with reconstructed change (dark: significant, p<0.05)',
                  fontsize=6.4)
    ax.tick_params(labelsize=6.4)
    ax.set_xlim(-max(downs) * 1.95, max(ups) * 1.22)
    for yy, u, d, su, sd in zip(y, ups, downs, sups, sdowns):
        ax.text(u + 200, yy, f'{u:,} ({su:,})', va='center', fontsize=5.0, color='#7B2C29')
        ax.text(-d - 150, yy, f'{d:,} ({sd:,})', va='center', ha='right', fontsize=5.0,
                color='#2C4E7B')
    handles = [Line2D([], [], marker='s', ls='', color=c, label=l, ms=5) for c, l in
               [('#F1C9C4', 'expanded'), ('#B03A2E', 'expanded sig.'),
                ('#CBDDF0', 'contracted'), ('#3B6EA5', 'contracted sig.')]]
    ax.legend(handles=handles, fontsize=5.0, loc='lower right', ncol=1,
              handletextpad=0.3)
    if own:
        ax.set_title('Family-size change per branch', fontsize=8.4, loc='left')
        save(fig, 'fig1b_cafe5_branch_bars')
    return ax


# ---------------------------------------------------------------- fig2
def load_scan():
    rows = list(csv.DictReader((SCAN / 'results_v2/clade_family_scan.tsv').open(), delimiter='\t'))
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
        p_up = float(r['p_cs'])
        base = float(r['cor']) + float(r['sol'])
        p_dn = binomtest(int(r['rac']), TOT['RAC'], min(1.0, base / (TOT['COR'] + TOT['SOL'])),
                         alternative='less').pvalue if base > 0 else 1.0
        cat, color = classify(r['symbol'], r['description'])
        out.append(dict(r, enr=enr, x=math.log2(enr), y=-math.log10(max(min(p_up, p_dn), 1e-300)),
                        p=min(p_up, p_dn), cat=cat, color=color))
    return out


def fig2(ax=None):
    rows = load_scan()
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for cat, color, _ in CATS + [('annotated (other)', '#BBBBBB', ''), ('unannotated', '#E6E6E6', '')]:
        xs = [r['x'] for r in rows if r['cat'] == cat]
        ys = [r['y'] for r in rows if r['cat'] == cat]
        if not xs:
            continue
        ax.scatter(xs, ys, s=3.2, c=color, alpha=0.6, lw=0, label=f'{cat} ({len(xs)})')
    ax.axvline(0, color='#888888', lw=0.6, ls='--')
    ax.axhline(-math.log10(0.05), color='#888888', lw=0.6, ls=':')
    ax.axvline(-math.log2(1.3), color='#BBBBBB', lw=0.5, ls='-.')
    ax.axvline(math.log2(1.3), color='#BBBBBB', lw=0.5, ls='-.')
    HIGHLIGHT = {'OG0000961': 'BLH1', 'OG0005156': 'ASY1', 'OG0005416': 'PDR2',
                 'OG0000864': 'CPK1', 'OG0000130': 'LNK1', 'OG0000222': 'RAPTOR2',
                 'OG0001188': 'Susy', 'OG0000664': 'enolase', 'OG0000348': 'MCM7',
                 'OG0005416 ': 'PDR2'}
    OFFSETS = {'OG0000256': (-6, 2), 'OG0000130': (-6, 2), 'OG0000222': (-6, 3),
               'OG0000664': (1, 6), 'OG0000864': (8, 5), 'OG0000348': (-9, 5),
               'OG0000961': (-9, -9), 'OG0005416': (4, 3), 'OG0001188': (-8, -7),
               'OG0005156': (6, -9)}
    for r in rows:
        if r['og'] in HIGHLIGHT:
            dx, dy = OFFSETS.get(r['og'], (3, 2))
            ax.annotate(HIGHLIGHT[r['og']], (r['x'], r['y']), xytext=(dx, dy),
                        textcoords='offset points',
                        ha='left' if dx > 0 else 'right', va='bottom',
                        fontsize=5.8, color='#7B2C29',
                        bbox=dict(fc='white', ec='none', alpha=0.6, pad=0.1),
                        arrowprops=dict(arrowstyle='-', lw=0.4, color='#999999',
                                        shrinkA=0, shrinkB=1))
    ax.set_xlabel('log2 enrichment (racemose vs diploid Prunus)', fontsize=7.2)
    ax.set_ylabel('\u2212log10 p', fontsize=7.2)
    ax.tick_params(labelsize=6.4)
    ax.legend(fontsize=5.0, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3,
              columnspacing=0.7, handletextpad=0.2, labelspacing=0.35)
    if own:
        n_all = len(rows)
        n_up = sum(1 for r in rows if r['x'] > 0 and r['p'] < 0.05)
        n_dn = sum(1 for r in rows if r['x'] < 0 and r['p'] < 0.05)
        n_cons = sum(1 for r in rows if float(r['p_cs']) < 0.05
                     and float(r['enr_cs']) > 1.3 and int(r['n_species_up']) >= 2
                     and r['te'] != 'True')
        n_te = sum(1 for r in rows if r['cat'] == 'TE / transposon')
        ax.text(0.015, 0.985,
                f'{n_all:,} families tested\n'
                f'{n_up:,} retained (p < 0.05)\n'
                f'{n_dn:,} fractionated (p < 0.05)\n'
                f'{n_cons} consistent\n'
                '(p < 0.05, >1.3x, >=2/3 spp., non-TE)\n'
                f'{n_te} TE-flagged (orange)\n'
                'dashed: 1x; dash-dot: 1.3x; dotted: p = 0.05',
                transform=ax.transAxes, fontsize=4.9, va='top', color='#555555',
                linespacing=1.45)
        ax.set_title('b  Genome-wide family scan by category', fontsize=8.4, loc='left')
        save(fig, 'fig2_scan_volcano')
    return ax


# ---------------------------------------------------------------- fig3
def fig3(ax=None):
    exp = list(csv.DictReader((SCAN / 'results_v2/racemose_expanded_consistent.tsv').open(),
                              delimiter='\t'))
    counts = {}
    with (B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups/'
              'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['Orthogroup']] = r
    tracks = ['pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B',
              'pcam', 'pav', 'pspe', 'pmal', 'par', 'pmu', 'pdu', 'psal', 'ppe', 'gtr', 'rch']
    glabels = ['PZI-A', 'PZI-B', 'PSP-A', 'PSP-B', 'PPA-A', 'PPA-B',
               'PCAM', 'PAV', 'PSPE', 'PMAL', 'PAR', 'PMU', 'PDU', 'PSAL', 'PPE', 'GTR', 'RCH']
    # pick top families, prioritising annotated non-TE candidates
    sel = [r for r in exp if r['symbol'] and float(r['rac']) >= 7 and 'transposon' not in
           (r['symbol'] + r['description']).lower() and 'retro' not in
           (r['symbol'] + r['description']).lower()]
    sel.sort(key=lambda r: -float(r['enr_cs']))
    top = sel[:12]
    M = np.zeros((len(top), len(tracks)))
    cats = []
    for i, r in enumerate(top):
        cr = counts[r['og']]
        for j, t in enumerate(tracks):
            M[i, j] = int(cr[t])
        cats.append(classify(r['symbol'], r['description']))
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
    im = ax.imshow(np.log1p(M), cmap='Reds', aspect='auto', vmin=0, vmax=np.log1p(M).max())
    ax.set_xticks(range(len(tracks)), glabels, rotation=90, fontsize=5.8)
    ax.set_yticks(range(len(top)), [f"{r['symbol']} ({r['og'][-5:]})" for r in top], fontsize=6.0)
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
    # category strip on the right (short codes + legend)
    code = {'TE / transposon': 'TE', 'floral / meristem / hormone': 'floral',
            'RNA / splicing': 'RNA', 'transport / vesicle': 'transport',
            'meiosis / DNA / chromatin': 'meiosis', 'metabolism / other': 'metab',
            'annotated (other)': 'other', 'unannotated': '-'}
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
    ax.legend(handles=handles, fontsize=5.4, loc='lower left', bbox_to_anchor=(0.62, 0.02),
              ncol=2, columnspacing=0.8, handletextpad=0.3)
    if own:
        ax.set_title('c  Consistent racemose-retained families (copy numbers; right = category)',
                     fontsize=8.4, loc='left', pad=16)
        save(fig, 'fig3_top_family_heatmap')
    return ax


# ---------------------------------------------------------------- fig4
def fig4(ax=None):
    rows = [r for r in csv.DictReader((SCAN / 'results_v2/go_enrichment_consistent.tsv').open(),
                                      delimiter='\t')]
    rows = [r for r in rows if float(r['q']) < 0.25][:8][::-1]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(4.6, 3.0))
    y = np.arange(len(rows))
    ax.scatter([-math.log10(max(float(r['q']), 1e-12)) for r in rows], y,
               s=[16 + 7 * int(r['expanded']) for r in rows], color='#C0504D', alpha=0.8, lw=0)
    for yy, r in zip(y, rows):
        ax.text(-math.log10(max(float(r['q']), 1e-12)) + 0.05, yy,
                f"{r['expanded']}/{r['bg']}", fontsize=5.6, va='center', color='#555555')
    ax.set_yticks(y, [r['go_description'][:34] for r in rows], fontsize=5.8)
    ax.set_xlabel('\u2212log10 FDR (BH)', fontsize=7.2)
    ax.tick_params(labelsize=6.4)
    ax.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    if own:
        ax.set_title('d  GO enrichment (racemose-retained)', fontsize=8.4, loc='left')
        save(fig, 'fig4_go_enrichment')
    return ax


# ---------------------------------------------------------------- fig5 (corrected model)
def fig5(ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(8.4, 3.9))
    ax.axis('off')
    ax.set_xlim(0, 13.2)
    ax.set_ylim(0, 6.6)

    def box(x, y, w, h, text, fc, ec='#555555', fs=5.9):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.07',
                                    fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs)

    def arrow(x1, y1, x2, y2, color='#555555'):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                     mutation_scale=8, color=color, lw=0.9))

    def chrom(cx, cy, h, color='#4F81BD'):
        ax.plot([cx, cx], [cy - h / 2, cy + h / 2], color=color, lw=2.0,
                solid_capstyle='round', zorder=4)
        ax.add_patch(Circle((cx, cy), 0.035, fc='white', ec=color, lw=0.6, zorder=5))

    def locus(x, y, n, label, color='#C0504D'):
        ax.plot([x, x + 0.52], [y, y], color='#999999', lw=1.0, zorder=3)
        for i in range(n):
            ax.add_patch(Circle((x + 0.15 + 0.22 * i, y), 0.045, fc=color, ec='none',
                                zorder=4))
        ax.text(x + 0.26, y - 0.24, label, ha='center', fontsize=4.9, color='#555555')

    # top note
    box(0.1, 6.02, 13.0, 0.5,
        'Shared WGD on the racemose stem only (after divergence from PCAM/PPE); '
        'PCAM and PPE are non-duplicated diploids.',
        '#FFFFFF', ec='#BBBBBB', fs=5.9)
    # ancestor
    box(0.15, 2.9, 2.05, 1.05, 'diploid Prunus ancestor\n(n = 8)', '#F2F2F2')
    icon_at(ax, 1.18, 4.55, 0.50, 'raceme')
    arrow(2.25, 3.85, 3.55, 4.85)
    arrow(2.25, 3.1, 3.55, 2.35)

    # racemose path
    box(3.6, 4.15, 3.2, 1.7,
        'racemose lineage\nPZI · PSP · PPA\n\nshared WGD (~55\u201342 Ma)\n'
        '11,747 expanded (956 significant)', '#FBE3D6')
    # 2n -> 4n schematic
    chrom(7.08, 4.92, 0.44); chrom(7.18, 4.92, 0.44)
    ax.text(7.13, 4.56, '2n', ha='center', fontsize=4.9, color='#555555')
    arrow(7.30, 4.92, 7.52, 4.92)
    for i, cx in enumerate([7.64, 7.74, 7.86, 7.96]):
        chrom(cx, 4.92, 0.44)
    ax.text(7.80, 4.56, '4n', ha='center', fontsize=4.9, color='#555555')
    locus(7.05, 4.22, 2, 'duplicates retained')
    icon_at(ax, 8.35, 4.9, 0.48, 'raceme')
    arrow(8.75, 4.9, 9.75, 4.9)
    box(9.85, 4.35, 3.2, 1.25, 'ancestral raceme\nprogramme retained\n(buffered duplicates)',
        '#DFF0D8')

    # diploid path
    box(3.6, 0.85, 3.2, 1.7,
        'corymbose & solitary lineages\nPCAM · PPE · \u2026\n\nno WGD\nno net copy-number gain',
        '#E8EEF7')
    chrom(7.08, 1.95, 0.44); chrom(7.18, 1.95, 0.44)
    ax.text(7.13, 1.59, '2n', ha='center', fontsize=4.9, color='#555555')
    locus(7.55, 1.95, 1, 'no net gain', color='#888888')
    icon_at(ax, 8.30, 1.9, 0.40, 'corymb')
    icon_at(ax, 9.12, 1.9, 0.40, 'solitary')
    arrow(9.45, 1.85, 9.75, 1.85)
    box(9.85, 1.1, 3.2, 1.25, 'derived architectures\n(copy-number stasis;\nregulatory divergence)',
        '#E8EEF7')

    ax.text(0.1, 0.28,
            'Control (CAFE5 without the WGD clade): diploid contractions disappear '
            '(rch/gtr/ppe \u2248 0 or positive) \u2014 the polarity is driven by '
            'racemose-side retention, not diploid loss.',
            fontsize=5.3, color='#555555')
    if own:
        ax.set_title('e  Working model: WGD-triggered retention of the ancestral racemose programme',
                     fontsize=8.4, loc='left')
        save(fig, 'fig5_evolution_model')
    return ax


def fig7(ax=None):
    """Shared vs clade-specific expansions: overlap of the three enrichment sets and their GO profiles."""
    counts = {}
    with (SCAN / 'results_v2/shared_specific_counts.tsv').open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['category']] = int(r['n_families'])
    rows = list(csv.DictReader((SCAN / 'results_v2/shared_specific_go.tsv').open(), delimiter='\t'))
    own = ax is None
    if own:
        fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.1),
                                 gridspec_kw={'width_ratios': [0.95, 1.3], 'wspace': 0.35})
    else:
        axes = [ax]
    axa = axes[0]
    axa.set_xlim(-1.55, 1.55)
    axa.set_ylim(-1.45, 1.25)
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
    axa.text(-1.05, 1.05, 'racemose\n(n=466)', color='#C0504D', fontsize=6.4, ha='center')
    axa.text(1.05, 1.05, 'corymbose\n(n=1,255)', color='#4F81BD', fontsize=6.4, ha='center')
    axa.text(0.0, -1.32, 'solitary (n=1,054)', color='#59A14F', fontsize=6.4, ha='center')
    if own:
        axa.set_title('a  Overlap of the enrichment sets', fontsize=8.0, loc='left')

    axb = axes[1]
    CAT_COL = {'R_only': '#C0504D', 'C_only': '#4F81BD', 'S_only': '#59A14F'}
    sel = []
    for g in CAT_COL:
        rs = sorted([r for r in rows if r['group'] == g], key=lambda r: float(r['q']))[:5]
        sel += rs[::-1]
    ys = np.arange(len(sel))
    for y, r in zip(ys, sel):
        col = CAT_COL[r['group']]
        x = -math.log10(max(float(r['q']), 1e-12))
        axb.plot([0, x], [y, y], color=col, lw=0.8, alpha=0.55)
        axb.scatter([x], [y], s=10 + 4 * int(r['k']), color=col, zorder=3)
    axb.set_yticks(ys, [r['go_description'][:38] for r in sel], fontsize=5.4)
    for tick, r in zip(axb.get_yticklabels(), sel):
        tick.set_color(CAT_COL[r['group']])
    axb.set_xlabel('\u2212log10 q (BH)', fontsize=7.0)
    axb.tick_params(axis='x', labelsize=6.4)
    axb.axvline(-math.log10(0.05), color='#888888', ls=':', lw=0.6)
    axb.set_xlim(0, max(-math.log10(max(float(r['q']), 1e-12)) for r in sel) * 1.12)
    handles = [Line2D([], [], marker='o', ls='', color=c, label=l, ms=4) for l, c in
               [('racemose-only', '#C0504D'), ('corymbose-only', '#4F81BD'),
                ('solitary-only', '#59A14F')]]
    axb.legend(handles=handles, fontsize=5.6, loc='lower right')
    if own:
        axb.set_title('b  GO profiles of the clade-specific sets', fontsize=8.0, loc='left')
        save(fig, 'fig7_shared_specific')
    return ax


def save(fig, name):
    for ext, kw in (('png', dict(dpi=600)), ('pdf', {}), ('svg', {})):
        target = FIG / f'{name}.{ext}'
        try:
            fig.savefig(target, bbox_inches='tight', **kw)
        except PermissionError:
            alt = FIG / f'{name}_v2.{ext}'
            fig.savefig(alt, bbox_inches='tight', **kw)
            print(f'WARNING: {target.name} is locked (open in a viewer?); wrote {alt.name}')
    fig.savefig(PREV / f'{name}_preview.png', dpi=150, bbox_inches='tight')
    print('saved', name)


# ---------------------------------------------------------------- fig0 integrated
# LEGACY: the current fig0_integrated is generated by plot_all17_figures.py
# (fig0_all17, all-17-track annotation). Kept for reference only; not called
# from __main__ so it cannot overwrite the current figure.
def fig0_legacy():
    fig = plt.figure(figsize=(11.5, 7.6))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.40,
                          left=0.06, right=0.985, top=0.93, bottom=0.07,
                          width_ratios=[1.0, 0.98, 0.95], height_ratios=[1.0, 1.05])

    # a framework: the dated tree (IQ-TREE + MCMCtree) with inflorescence types
    axa = fig.add_subplot(gs[0, 0])
    root, H, ids = load_tree()
    draw_chronogram(axa, root, H, tags={'<11>': ('WGD', '#B03A2E')},
                    xpad=46, fs_label=5.4, show_axis=False)
    for node in iter_nodes(root):
        nm = node['name']
        if nm in INFL_TYPE:
            icon_at(axa, H + 35, node['y'], 0.22, INFL_TYPE[nm])
    axa.set_title('a  Framework: dated tree, inflorescence types and the WGD',
                  fontsize=8.4, loc='left')

    axb = fig.add_subplot(gs[0, 1])
    fig1(axb, compact=True)
    axb.set_title('b  CAFE5 family-size change along the tree', fontsize=8.4, loc='left')

    axc = fig.add_subplot(gs[0, 2])
    fig2(axc)
    axc.set_title('c  Genome-wide family scan (by category)', fontsize=8.4, loc='left')

    axd = fig.add_subplot(gs[1, 0:2])
    fig3(axd)
    axd.set_title('d  Consistent racemose-retained families (right: functional category)',
                  fontsize=8.4, loc='left', pad=16)

    axe = fig.add_subplot(gs[1, 2])
    fig4(axe)
    axe.set_title('e  GO enrichment of retained families', fontsize=8.4, loc='left')

    save(fig, 'fig0_integrated')


if __name__ == '__main__':
    fig1()
    fig1_bars()
    fig2()
    fig3()
    fig4()
    fig5()
    fig7()
    # fig0 is generated by plot_all17_figures.py::fig0_all17 (all-17 annotation)
