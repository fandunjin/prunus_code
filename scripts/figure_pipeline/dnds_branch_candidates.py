#!/usr/bin/env python3
"""dnds_branch_candidates.py - branch-model dN/dS tests for inflorescence candidate
subfamilies across six species (pzi-A1, psp-A1, ppa-A, PCAM, PPE, GTR).

For each target subfamily the six-species family tree (inflorescence pipeline) is
used to extract one-sequence-per-species ortholog sets (preferred A1 haplotype for
the racemose species). Each set is aligned (mafft protein -> pal2nal codons) and
tested with PAML codeml branch models:
    null    : single omega
    branch  : racemose clade (stem + ppa/psp/pzi) as foreground omega
LRT = 2 dlnL, df = 1.

Outputs (analysis/dnds_branch_20260922/):
    results/branch_dnds_results.tsv     per ortholog set
    results/branch_dnds_summary.tsv     per subfamily
    runs/<set>/                         alignments, trees, codeml output
"""
import csv
import hashlib
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from scipy.stats import chi2

B = Path('${PROJECT_ROOT}')
PIPE = B / 'analysis/inflorescence_family_pipeline_crossspecies6_20260913'
HAPL = B / '0903/haploid_repartition/inputs'
PREP = B / '0903/karyotype/prepared'
OUT = B / 'analysis/dnds_branch_20260922'
RES = OUT / 'results'
RUNS = OUT / 'runs'
RES.mkdir(parents=True, exist_ok=True)
RUNS.mkdir(parents=True, exist_ok=True)

MAFFT = '${CONDA_ROOT}/bin/mafft'
PAL2NAL = '${CONDA_ROOT}/bin/pal2nal.pl'
CODEML = '${CONDA_ROOT}/bin/codeml'

SPECIES_FILES = {
    'pzi': (HAPL / 'pzi/pep/pzi_hapA1.pep.fa', HAPL / 'pzi/cds/pzi_hapA1.cds.fa'),
    'psp': (HAPL / 'psp/pep/psp_hapA1.pep.fa', HAPL / 'psp/cds/psp_hapA1.cds.fa'),
    'ppa': (HAPL / 'ppa/pep/ppa_hapA.pep.fa', HAPL / 'ppa/cds/ppa_hapA.cds.fa'),
    'pcam': (B / 'analysis/canonical/pcam/primary.pep.fa', B / 'analysis/canonical/pcam/primary.cds.fa'),
    'ppe': (B / 'analysis/canonical/ppe/primary.pep.fa', B / 'analysis/canonical/ppe/primary.cds.fa'),
    'gtr': (PREP / 'gtr/gtr.pep.fa', PREP / 'gtr/gtr.cds.fa'),
}

TARGETS = [  # (family, subfamily_group, candidate?)
    ('MADS', 'SVP/DAM', True),
    ('MADS', 'AP1/CAL', False),
    ('MADS', 'AGL6', False),
    ('MADS', 'SEP', False),
    ('MADS', 'SOC1', False),
    ('MADS', 'AG/SHP/STK', False),
    ('KNOX', 'KNAT1/2/6', True),
    ('KNOX', 'STM', True),
    ('YABBY', 'YABBY1/2/3', True),
    ('YABBY', 'YABBY4/5', False),
    ('SPL', 'SPL9/15', True),
    ('SPL', 'SPL3/4/5', False),
    ('SPL', 'SPL-other', False),
]

TREE6 = '(gtr:0.3, ((ppa:0.1 #1, (psp:0.1 #1, pzi:0.1 #1):0.05 #1):0.05 #1, (pcam:0.1, ppe:0.1):0.1):0.2);'
TREE6_NULL = '(gtr:0.3, ((ppa:0.1, (psp:0.1, pzi:0.1):0.05):0.05, (pcam:0.1, ppe:0.1):0.1):0.2);'


def read_fasta(path):
    d, name, buf = {}, None, []
    for line in Path(path).read_text(errors='ignore').splitlines():
        if line.startswith('>'):
            if name:
                d[name] = ''.join(buf)
            name = line[1:].split()[0]
            buf = []
        elif name:
            buf.append(line.strip())
    if name:
        d[name] = ''.join(buf)
    return d


def md5(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def load_species():
    """species -> {'by_id': id->cds, 'by_md5': md5(pep)->(id, cds)}"""
    out = {}
    for sp, (pep_f, cds_f) in SPECIES_FILES.items():
        pep, cds = read_fasta(pep_f), read_fasta(cds_f)

        def norm(x):
            x = re.sub(r'-mRNA-\d+$', '', x)
            x = re.sub(r'\.p$', '', x)
            return x

        norm_cds = {}
        for cid, cseq in cds.items():
            norm_cds[norm(cid)] = (cid, cseq)
            norm_cds[cid] = (cid, cseq)
        by_md5 = {}
        for pid, seq in pep.items():
            key = md5(seq.rstrip('*'))
            if norm(pid) in norm_cds:
                by_md5[key] = norm_cds[norm(pid)]
            elif pid in norm_cds:
                by_md5[key] = norm_cds[pid]
        out[sp] = {'by_md5': by_md5, 'by_id': norm_cds}
    return out


def get_cds(sp_data, protein):
    s = protein.rstrip('*')
    hit = sp_data['by_md5'].get(md5(s))
    if hit:
        return hit
    return None


# ---------------------------------------------------------------- tree tools
def parse_newick(text):
    s = text.strip()
    pos = 0

    def parse_node():
        nonlocal pos
        node = {'children': [], 'name': '', 'length': 0.0}
        if s[pos] == '(':
            pos += 1
            while True:
                node['children'].append(parse_node())
                if s[pos] == ',':
                    pos += 1
                    continue
                if s[pos] == ')':
                    pos += 1
                    break
            m = re.match(r'([^:,()\s]+)', s[pos:])
            if m and s[pos:pos + 1] not in ':,)':
                node['name'] = m.group(1)
                pos += len(node['name'])
        else:
            m = re.match(r'([^:,()\s]+)', s[pos:])
            if m:
                node['name'] = m.group(1)
                pos += len(node['name'])
        if pos < len(s) and s[pos] == ':':
            pos += 1
            m = re.match(r'([\d.eE+-]+)', s[pos:])
            if m:
                node['length'] = float(m.group(1))
                pos += len(m.group(1))
        return node

    return parse_node()


def leaves_of(n):
    if not n['children']:
        return [n]
    out = []
    for c in n['children']:
        out += leaves_of(c)
    return out


def all_nodes(root):
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack += n['children']
    return out


def spe_of(name):
    return name.split('_', 1)[0]


def find_ortholog_sets(root, pref_tips, prot, max_sets=4):
    """pref_tips: species -> {normalized tip name: full tip}. Returns list of
    dicts species -> tip (one tip per species)."""
    tip_by_norm = {}
    for sp, tips in pref_tips.items():
        for t in tips:
            tip_by_norm[t] = sp
    sets = []
    for node in all_nodes(root):
        lv = []
        for l in leaves_of(node):
            nm = re.sub(r'-mRNA-\d+$', '', l['name'])
            if nm in tip_by_norm:
                lv.append(nm)
        if not lv:
            continue
        by_sp = defaultdict(list)
        for t in lv:
            by_sp[tip_by_norm[t]].append(t)
        if len(by_sp) < 6:
            continue
        chosen = {}
        for sp, tips in by_sp.items():
            chosen[sp] = max(tips, key=lambda t: len(prot.get(t, '')))
        # avoid sets that keep growing (take smallest clade containing all six)
        sets.append((len(lv), chosen))
    # keep smallest clades (most specific ortholog sets), dedupe by membership
    sets.sort(key=lambda x: x[0])
    out, seen = [], set()
    for _, chosen in sets:
        key = tuple(sorted(chosen.values()))
        if key in seen:
            continue
        seen.add(key)
        out.append(chosen)
        if len(out) >= max_sets:
            break
    return out


# ---------------------------------------------------------------- codeml
PARAMS = """noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = {model}
NSsites = 0
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = .4
fix_alpha = 1
alpha = 0.
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
fix_blength = 2
method = 0
"""


def run_codeml(run_dir, aln_phy, tree_text, model, tag):
    (run_dir / f'tree_{tag}.nwk').write_text(tree_text + '\n')
    (run_dir / f'codeml_{tag}.ctl').write_text(
        f'seqfile = {aln_phy}\ntreefile = {run_dir}/tree_{tag}.nwk\n'
        f'outfile = {run_dir}/mlc_{tag}\n' + PARAMS.format(model=model))
    r = subprocess.run([CODEML, f'codeml_{tag}.ctl'], cwd=run_dir,
                       capture_output=True, text=True)
    out = run_dir / f'mlc_{tag}'
    if not out.exists():
        return None
    txt = out.read_text(errors='ignore')
    m = re.search(r'lnL\(ntime.*?\):\s*(-?[\d.]+)', txt)
    if not m:
        return None
    lnL = float(m.group(1))
    om = re.findall(r'omega \(dN/dS\) =\s*([\d.]+)', txt)
    omb = re.findall(r'w \(dN/dS\) for branches:\s*([\d.]+)\s+([\d.]+)', txt)
    omegas = [float(x) for x in om]
    if omb:
        omegas = [float(omb[0][0]), float(omb[0][1])]
    if model == 0:
        return {'lnL': lnL, 'omega_bg': omegas[0] if omegas else None, 'omega_fg': None}
    return {'lnL': lnL,
            'omega_bg': omegas[0] if omegas else None,
            'omega_fg': omegas[1] if len(omegas) > 1 else None}


def main():
    sp_data = load_species()
    cand = defaultdict(lambda: defaultdict(list))
    with (PIPE / 'candidates_crossspecies.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            key = (r['family'], r['subfamily_group'])
            cand[key][(r['species'], r['phase'])].append(r['gene_key'])

    results = []
    for family, sub, is_cand in TARGETS:
        tree_file = PIPE / 'trees' / f'{family}.treefile'
        if not tree_file.exists():
            print('missing tree', family)
            continue
        fam_fa = read_fasta(PIPE / 'sequences' / f'{family}.fa')
        prot = {re.sub(r'-mRNA-\d+$', '', k): v for k, v in fam_fa.items()}
        root = parse_newick(tree_file.read_text(errors='ignore'))
        # preferred tips
        pref_species = {'pzi': 'A1', 'psp': 'A1', 'ppa': 'A', 'pcam': 'S', 'ppe': 'S', 'gtr': 'S'}
        pref_tips = defaultdict(dict)
        for (sp, ph), keys in cand[(family, sub)].items():
            if pref_species.get(sp) != ph:
                continue
            for gk in keys:
                pref_tips[sp][re.sub(r'-mRNA-\d+$', '', gk)] = gk
        orth_sets = find_ortholog_sets(root, pref_tips, prot)
        print(f'{family} | {sub}: {len(orth_sets)} ortholog set(s)', flush=True)
        for si, chosen in enumerate(orth_sets):
            tag = f'{family}_{sub.replace("/", "-")}_{si}'
            run_dir = RUNS / tag
            run_dir.mkdir(parents=True, exist_ok=True)
            pep_lines, cds_lines, ok = [], [], True
            for sp in ['pzi', 'psp', 'ppa', 'pcam', 'ppe', 'gtr']:
                gk = chosen.get(sp)
                if not gk:
                    ok = False
                    break
                protein = prot.get(gk) or fam_fa.get(gk)
                if not protein:
                    ok = False
                    break
                hit = get_cds(sp_data[sp], protein)
                if not hit:
                    ok = False
                    print(f'  {tag}: no CDS for {sp} {gk}')
                    break
                cid, cseq = hit
                pep_lines.append(f'>{sp}\n{protein}')
                cds_lines.append(f'>{sp}\n{cseq}')
            if not ok:
                continue
            (run_dir / 'pep.fa').write_text('\n'.join(pep_lines) + '\n')
            (run_dir / 'cds.fa').write_text('\n'.join(cds_lines) + '\n')
            with (run_dir / 'aln.fa').open('w') as fh:
                subprocess.run([MAFFT, '--auto', '--quiet', str(run_dir / 'pep.fa')],
                               stdout=fh, check=True)
            r = subprocess.run([PAL2NAL, str(run_dir / 'aln.fa'), str(run_dir / 'cds.fa'),
                                '-output', 'fasta', '-nogap'],
                               capture_output=True, text=True)
            if r.returncode != 0 or len(r.stdout) < 200:
                print(f'  {tag}: pal2nal failed')
                continue
            (run_dir / 'codon.fa').write_text(r.stdout)
            # to phylip
            seqs = read_fasta(run_dir / 'codon.fa')
            seqs = {k: v for k, v in seqs.items() if k in ('pzi', 'psp', 'ppa', 'pcam', 'ppe', 'gtr')}
            if len(seqs) < 6:
                continue
            L = len(next(iter(seqs.values())))
            with (run_dir / 'aln.phy').open('w') as fh:
                fh.write(f'{len(seqs)} {L}\n')
                for k, v in seqs.items():
                    fh.write(f'{k:<12s}{v}\n')
            null = run_codeml(run_dir, run_dir / 'aln.phy', TREE6_NULL, 0, 'null')
            alt = run_codeml(run_dir, run_dir / 'aln.phy', TREE6, 2, 'branch')
            if not null or not alt:
                print(f'  {tag}: codeml failed')
                continue
            dlnL = 2 * (alt['lnL'] - null['lnL'])
            p = float(chi2.sf(max(dlnL, 0), 1))
            results.append(dict(family=family, subfamily=sub, candidate=is_cand,
                                set=tag, n_codons=L,
                                omega_bg=null['omega_bg'],
                                omega_rac=alt['omega_fg'],
                                omega_bg_branch=alt['omega_bg'],
                                lnL_null=null['lnL'], lnL_branch=alt['lnL'],
                                LRT=dlnL, p=p))
            print(f'  {tag}: L={L} w_bg={null["omega_bg"]} w_rac={alt["omega_fg"]} p={p:.3g}',
                  flush=True)

    cols = list(results[0].keys()) if results else []
    with (RES / 'branch_dnds_results.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow(r)
    # summary per subfamily
    by = defaultdict(list)
    for r in results:
        by[(r['family'], r['subfamily'])].append(r)
    with (RES / 'branch_dnds_summary.tsv').open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['family', 'subfamily', 'candidate', 'n_sets', 'omega_rac_mean',
                    'omega_bg_mean', 'ratio', 'min_p'])
        for key, rs in sorted(by.items()):
            wr = [r['omega_rac'] for r in rs if r['omega_rac'] is not None]
            wb = [r['omega_bg'] for r in rs if r['omega_bg'] is not None]
            w.writerow([key[0], key[1], rs[0]['candidate'], len(rs),
                        f'{sum(wr)/len(wr):.3f}' if wr else '',
                        f'{sum(wb)/len(wb):.3f}' if wb else '',
                        f'{(sum(wr)/len(wr))/(sum(wb)/len(wb)):.2f}' if wr and wb else '',
                        f'{min(r["p"] for r in rs):.3g}'])
    print('written', RES)


if __name__ == '__main__':
    main()
