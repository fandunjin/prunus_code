#!/usr/bin/env python3
"""dnds_background_controls.py - genome-wide background for the racemose branch
dN/dS test. Random single-copy orthogroups (present in all six tracks) are
processed exactly like the candidate sets: six-species codon alignment +
codeml branch model with the racemose clade as foreground.

Outputs (results/):
    branch_dnds_background.tsv
"""
import csv
import hashlib
import random
import re
import subprocess
from pathlib import Path

from scipy.stats import chi2

B = Path('${PROJECT_ROOT}')
FROZEN = B / '0903/analysisphylogenetics/01_frozen_inputs/main/pep'
ORTHO = B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups'
HAPL = B / '0903/haploid_repartition/inputs'
CAN = B / 'analysis/canonical'
PREP = B / '0903/karyotype/prepared'
OUT = B / 'analysis/dnds_branch_20260922'
RES = OUT / 'results'
RUNS = OUT / 'runs_bg_B'
RUNS.mkdir(parents=True, exist_ok=True)

MAFFT = '${CONDA_ROOT}/bin/mafft'
PAL2NAL = '${CONDA_ROOT}/bin/pal2nal.pl'
CODEML = '${CONDA_ROOT}/bin/codeml'

TRACKS = {'pzi': 'pzi_B.fa', 'psp': 'psp_B.fa', 'ppa': 'ppa_B.fa',
          'pcam': 'pcam.fa', 'ppe': 'ppe.fa', 'gtr': 'gtr.fa'}
SPECIES_FILES = {
    'pzi': (HAPL / 'pzi/pep/pzi_hapB1.pep.fa', HAPL / 'pzi/cds/pzi_hapB1.cds.fa'),
    'psp': (HAPL / 'psp/pep/psp_hapB1.pep.fa', HAPL / 'psp/cds/psp_hapB1.cds.fa'),
    'ppa': (HAPL / 'ppa/pep/ppa_hapB.pep.fa', HAPL / 'ppa/cds/ppa_hapB.cds.fa'),
    'pcam': (CAN / 'pcam/primary.pep.fa', CAN / 'pcam/primary.cds.fa'),
    'ppe': (CAN / 'ppe/primary.pep.fa', CAN / 'ppe/primary.cds.fa'),
    'gtr': (PREP / 'gtr/gtr.pep.fa', PREP / 'gtr/gtr.cds.fa'),
}

TREE6 = '(gtr:0.3, ((ppa:0.1 #1, (psp:0.1 #1, pzi:0.1 #1):0.05 #1):0.05 #1, (pcam:0.1, ppe:0.1):0.1):0.2);'
TREE6_NULL = '(gtr:0.3, ((ppa:0.1, (psp:0.1, pzi:0.1):0.05):0.05, (pcam:0.1, ppe:0.1):0.1):0.2);'

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


def md5(s):
    return hashlib.md5(s.rstrip('*').encode()).hexdigest()


def load_md5_cds():
    out = {}
    for sp, (pep_f, cds_f) in SPECIES_FILES.items():
        pep, cds = read_fasta(pep_f), read_fasta(cds_f)

        def norm(x):
            x = re.sub(r'-mRNA-\d+$', '', x)
            x = re.sub(r'\.p$', '', x)
            return x

        norm_cds = {}
        for cid, cseq in cds.items():
            norm_cds[norm(cid)] = cseq
            norm_cds[cid] = cseq
        by = {}
        for pid, seq in pep.items():
            cseq = norm_cds.get(norm(pid)) or norm_cds.get(pid)
            if cseq:
                by[md5(seq)] = cseq
        out[sp] = by
        print(sp, 'md5->cds entries:', len(by))
    return out


def run_codeml(run_dir, aln_phy, tree_text, model, tag):
    (run_dir / f'tree_{tag}.nwk').write_text(tree_text + '\n')
    (run_dir / f'codeml_{tag}.ctl').write_text(
        f'seqfile = {aln_phy}\ntreefile = {run_dir}/tree_{tag}.nwk\n'
        f'outfile = {run_dir}/mlc_{tag}\n' + PARAMS.format(model=model))
    subprocess.run([CODEML, f'codeml_{tag}.ctl'], cwd=run_dir, capture_output=True, text=True)
    out = run_dir / f'mlc_{tag}'
    if not out.exists():
        return None
    txt = out.read_text(errors='ignore')
    m = re.search(r'lnL\(ntime.*?\):\s*(-?[\d.]+)', txt)
    if not m:
        return None
    lnL = float(m.group(1))
    if model == 0:
        om = re.findall(r'omega \(dN/dS\) =\s*([\d.]+)', txt)
        return {'lnL': lnL, 'omega_bg': float(om[0]) if om else None, 'omega_fg': None}
    omb = re.findall(r'w \(dN/dS\) for branches:\s*([\d.]+)\s+([\d.]+)', txt)
    if not omb:
        return None
    return {'lnL': lnL, 'omega_bg': float(omb[0][0]), 'omega_fg': float(omb[0][1])}


def main():
    md5cds = load_md5_cds()
    # single-copy orthogroups present in the six chosen tracks
    counts = {}
    with (ORTHO / 'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            counts[r['Orthogroup']] = r
    tracks = ['pzi_B', 'psp_B', 'ppa_B', 'pcam', 'ppe', 'gtr']
    candidates = []
    with (ORTHO / 'Orthogroups_SingleCopyOrthologues.txt').open(errors='ignore') as fh:
        for line in fh:
            og = line.strip()
            r = counts.get(og)
            if not r:
                continue
            if all(int(r[t]) == 1 for t in tracks):
                candidates.append(og)
    print('single-copy OGs in the six tracks:', len(candidates))
    random.seed(20260922)
    sample = random.sample(candidates, min(40, len(candidates)))

    genes = {}
    with (ORTHO / 'Orthogroups.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            genes[r['Orthogroup']] = r

    seqs = {sp: read_fasta(FROZEN / fn) for sp, fn in TRACKS.items()}

    results = []
    for og in sample:
        row = genes.get(og)
        if not row:
            continue
        run_dir = RUNS / og
        run_dir.mkdir(parents=True, exist_ok=True)
        pep_lines, cds_lines, ok = [], [], True
        for sp, track in [('pzi', 'pzi_B'), ('psp', 'psp_B'), ('ppa', 'ppa_B'),
                          ('pcam', 'pcam'), ('ppe', 'ppe'), ('gtr', 'gtr')]:
            gid = (row.get(track, '') or '').split(',')[0].strip()
            protein = seqs[sp].get(gid)
            if not protein:
                ok = False
                break
            cseq = md5cds[sp].get(md5(protein))
            if not cseq:
                ok = False
                break
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
                            '-output', 'fasta', '-nogap'], capture_output=True, text=True)
        if r.returncode != 0 or len(r.stdout) < 200:
            continue
        (run_dir / 'codon.fa').write_text(r.stdout)
        seqd = read_fasta(run_dir / 'codon.fa')
        seqd = {k: v for k, v in seqd.items() if k in ('pzi', 'psp', 'ppa', 'pcam', 'ppe', 'gtr')}
        if len(seqd) < 6:
            continue
        L = len(next(iter(seqd.values())))
        with (run_dir / 'aln.phy').open('w') as fh:
            fh.write(f'{len(seqd)} {L}\n')
            for k, v in seqd.items():
                fh.write(f'{k:<12s}{v}\n')
        null = run_codeml(run_dir, run_dir / 'aln.phy', TREE6_NULL, 0, 'null')
        alt = run_codeml(run_dir, run_dir / 'aln.phy', TREE6, 2, 'branch')
        if not null or not alt or alt['omega_fg'] is None:
            continue
        dlnL = 2 * (alt['lnL'] - null['lnL'])
        p = float(chi2.sf(max(dlnL, 0), 1))
        results.append(dict(og=og, n_codons=L, omega_bg=null['omega_bg'],
                            omega_rac=alt['omega_fg'], omega_bg_branch=alt['omega_bg'],
                            LRT=dlnL, p=p))
        print(f'{og}: L={L} w_bg={null["omega_bg"]} w_rac={alt["omega_fg"]} p={p:.3g}', flush=True)

    with (RES / 'branch_dnds_background_B.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    ps = [r['p'] for r in results]
    sig = sum(1 for p in ps if p < 0.05)
    print(f'background: n={len(ps)}  p<0.05: {sig} ({sig/len(ps):.1%})')
    print('written', RES / 'branch_dnds_background_B.tsv')


if __name__ == '__main__':
    main()
