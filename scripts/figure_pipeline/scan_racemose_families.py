#!/usr/bin/env python3
"""scan_racemose_families.py - genome-wide gene-family scan for racemose-specific
retention/expansion relative to the corymbose and solitary clades and to the
non-Prunus outgroups.

Groups (17 OrthoFinder tracks)
  racemose     : ppa_A, ppa_B, psp_A, psp_B, pzi_A, pzi_B   (shared WGD, x2 expected)
  corymbose    : pcam, pav, pspe, pmal                      (Cerasus clade)
  solitary     : par, pmu, pdu, psal, ppe                   (Prunus/Amygdalus clade)
  outgroup     : gtr, rch

Test: the racemose tracks carry two copies per ancestral locus (WGD); for each
family the expected racemose count under one-to-one retention is
    2 * p_diploid * N_racemose_genes,  p_diploid = family copies / diploid genes
tested with a binomial (one-sided, greater and less). Deviations flag
racemose-specific expansions (retention above the WGD expectation) and
depletions (fractionation).

Annotation: racemose-track member genes are mapped back to the original gene IDs
(frozen id_maps) and annotated with the Swiss-Prot transfer tables (pzi/psp/ppa).

Outputs (results/)
  racemose_family_scan.tsv        all tested families with counts, enrichment, p
  racemose_expanded_top.tsv       significantly expanded families (ranked)
  racemose_depleted_top.tsv       significantly depleted families (ranked)
  racemose_unique_families.tsv    families with no diploid-Prunus representative
  go_enrichment_expanded.tsv      GO terms over-represented among expanded families
"""
import csv
import math
import re
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import binomtest, fisher_exact


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        prev = min(prev, pvals[idx] * n / i)
        q[idx] = prev
    return q

B = Path('${PROJECT_ROOT}')
ORTHO = B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups'
FROZEN = B / '0903/analysisphylogenetics/01_frozen_inputs/main/id_maps'
ANN = B / '0903/annotation/functional'
OUT = B / 'analysis/racemose_inflorescence_scan_20260922/results'
OUT.mkdir(parents=True, exist_ok=True)

RAC = ['ppa_A', 'ppa_B', 'psp_A', 'psp_B', 'pzi_A', 'pzi_B']
COR = ['pcam', 'pav', 'pspe', 'pmal']
SOL = ['par', 'pmu', 'pdu', 'psal', 'ppe']
OUTG = ['gtr', 'rch']

# ---------------------------------------------------------------- counts
counts = {}          # og -> {track: n}
totals = Counter()
with (ORTHO / 'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    for r in rd:
        og = r['Orthogroup']
        d = {t: int(r[t]) for t in RAC + COR + SOL + OUTG}
        counts[og] = d
        for t, v in d.items():
            totals[t] += v

TOT = {
    'RAC': sum(totals[t] for t in RAC),
    'COR': sum(totals[t] for t in COR),
    'SOL': sum(totals[t] for t in SOL),
    'OUT': sum(totals[t] for t in OUTG),
    'COR_SOL': sum(totals[t] for t in COR + SOL),
}
print('gene totals:', dict(TOT), 'families:', len(counts))


def summ(d, tracks):
    return sum(d[t] for t in tracks)


def test_family(obs, n_total, base_count, base_total):
    """Binomial test of the racemose copy number against the one-to-one
    (per-track density) expectation from a diploid baseline. The racemose tracks
    already contain both A and B homeologs, so the uniform WGD doubling is the
    null; deviations flag family-specific retention (>1) or fractionation (<1)
    on the racemose lineage."""
    if base_total == 0 or n_total == 0:
        return None
    p = min(1.0, base_count / base_total)
    exp = p * n_total
    pg = binomtest(obs, n_total, p, alternative='greater').pvalue
    pl = binomtest(obs, n_total, p, alternative='less').pvalue
    return exp, pg, pl


# ---------------------------------------------------------------- annotation
def load_id_map(track):
    m = {}
    f = FROZEN / f'{track}.id_map.tsv'
    if not f.exists():
        return m
    with f.open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            m[r['stable_id']] = r['original_id']
    return m


IDMAP = {t: load_id_map(t) for t in RAC}

DESC = {}          # species -> protein_id -> dict(description, symbols, go)
for sp in ['pzi', 'psp', 'ppa']:
    d = {}
    f = ANN / sp / 'swissprot_functional_annotation.tsv'
    if f.exists():
        with f.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                d[r['protein_id']] = {
                    'desc': (r.get('protein_description') or '').split(';')[0][:140],
                    'symbol': (r.get('gene_names') or '').split()[0] if r.get('gene_names') else '',
                    'go': r.get('go_terms') or '',
                }
    DESC[sp] = d

OG_GENES = {}
with (ORTHO / 'Orthogroups.tsv').open(errors='ignore') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        OG_GENES[r['Orthogroup']] = r


def annotate_og(og):
    row = OG_GENES.get(og, {})
    symbols, descs, gos = Counter(), [], set()
    for tr in RAC:
        sp = tr.split('_')[0]
        for gid in [x.strip() for x in (row.get(tr, '') or '').split(',') if x.strip()]:
            orig = IDMAP.get(tr, {}).get(gid, gid)
            a = DESC.get(sp, {}).get(orig)
            if not a:
                continue
            if a['symbol']:
                symbols[a['symbol']] += 1
            if a['desc'] and len(descs) < 3:
                descs.append(a['desc'])
            for term in re.split(r'[;,]', a['go']):
                term = term.strip()
                if term.startswith('GO:'):
                    gos.add(term)
    return (symbols.most_common(1)[0][0] if symbols else '',
            ' | '.join(dict.fromkeys(descs))[:160], gos)


# ---------------------------------------------------------------- scan
rows = []
for og, d in counts.items():
    rac = summ(d, RAC)
    cor = summ(d, COR)
    sol = summ(d, SOL)
    out = summ(d, OUTG)
    present = sum(1 for t in RAC + COR + SOL + OUTG if d[t] > 0)
    if present < 2 or rac + cor + sol + out > 100:
        continue
    if cor + sol < 3:          # need an informed diploid baseline
        continue
    res = {}
    for label, cnt, tot in [('cor', cor, TOT['COR']), ('sol', sol, TOT['SOL']),
                            ('cor_sol', cor + sol, TOT['COR_SOL']), ('out', out, TOT['OUT'])]:
        r = test_family(rac, TOT['RAC'], cnt, tot)
        if r is None:
            res[label] = (float('nan'), 1.0, 1.0)
        else:
            exp, pg, pl = r
            res[label] = (rac / exp if exp > 0 else float('inf'), pg, pl)
    rows.append((og, rac, cor, sol, out, res))

# annotate only families that pass a coarse screen (speed)
keep = [r for r in rows
        if (r[5]['cor_sol'][1] < 0.05 and r[5]['cor_sol'][0] > 1.3)
        or (r[5]['cor_sol'][2] < 0.05 and r[5]['cor_sol'][0] < 0.7)
        or (r[1] > 0 and r[2] + r[3] == 0)]
print('families passing screen for annotation:', len(keep))

annot = {r[0]: annotate_og(r[0]) for r in keep}

with (OUT / 'racemose_family_scan.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['og', 'rac_copies', 'cor_copies', 'sol_copies', 'out_copies',
                'enr_vs_cor', 'p_cor_up', 'p_cor_down',
                'enr_vs_sol', 'p_sol_up', 'p_sol_down',
                'enr_vs_cor_sol', 'p_up', 'p_down',
                'enr_vs_out', 'p_out_up', 'p_out_down',
                'symbol', 'description'])
    for og, rac, cor, sol, out, res in rows:
        e_cs, p_cs_up, p_cs_dn = res['cor_sol']
        sym, desc, gos = annot.get(og, ('', '', set()))
        w.writerow([og, rac, cor, sol, out,
                    f'{res["cor"][0]:.3f}', f'{res["cor"][1]:.3g}', f'{res["cor"][2]:.3g}',
                    f'{res["sol"][0]:.3f}', f'{res["sol"][1]:.3g}', f'{res["sol"][2]:.3g}',
                    f'{e_cs:.3f}', f'{p_cs_up:.3g}', f'{p_cs_dn:.3g}',
                    f'{res["out"][0]:.3f}', f'{res["out"][1]:.3g}', f'{res["out"][2]:.3g}',
                    sym, desc])

# top expanded / depleted by the pooled diploid-Prunus test
def collect(pred, reverse=False):
    out = []
    for og, rac, cor, sol, out_c, res in rows:
        e, pu, pd = res['cor_sol']
        if pred(e, pu, pd):
            sym, desc, gos = annot.get(og, ('', '', set()))
            out.append((og, rac, cor, sol, out_c, e, pu, pd, sym, desc))
    out.sort(key=lambda x: -x[5] if reverse is False else x[5])
    return out

expanded = collect(lambda e, pu, pd: pu < 0.05 and e > 1.3)
depleted = collect(lambda e, pu, pd: pd < 0.05 and e < 0.7)
unique = [(og, r[1], r[2], r[3], r[4]) for og, r in
          ((og, (og, rac, cor, sol, out)) for og, rac, cor, sol, out, _ in rows)
          if r[2] + r[3] == 0 and r[1] > 0]

with (OUT / 'racemose_expanded_top.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['og', 'rac_copies', 'cor_copies', 'sol_copies', 'out_copies',
                'enr_vs_cor_sol', 'p_up', 'p_down', 'symbol', 'description'])
    for r in expanded:
        w.writerow(r)
with (OUT / 'racemose_depleted_top.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['og', 'rac_copies', 'cor_copies', 'sol_copies', 'out_copies',
                'enr_vs_cor_sol', 'p_up', 'p_down', 'symbol', 'description'])
    for r in depleted:
        w.writerow(r)
with (OUT / 'racemose_unique_families.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['og', 'rac_copies', 'cor_copies', 'sol_copies', 'out_copies'])
    for r in unique:
        w.writerow(r)

print(f'expanded: {len(expanded)}  depleted: {len(depleted)}  racemose-unique: {len(unique)}')

# ---------------------------------------------------------------- GO enrichment
# family -> GO set (annotated subset only)
fam_go = {}
for og, (sym, desc, gos) in annot.items():
    if gos:
        fam_go[og] = gos
exp_set = {r[0] for r in expanded if r[0] in fam_go}
bg_set = set(fam_go)
print(f'GO-annotated families: {len(bg_set)}; expanded with GO: {len(exp_set)}')

go_count = defaultdict(int)
for og in exp_set:
    for g in fam_go[og]:
        go_count[g] += 1
results = []
for term, k in go_count.items():
    if k < 3:
        continue
    bg_k = sum(1 for og in bg_set if term in fam_go[og])
    table = [[k, len(exp_set) - k], [bg_k - k, len(bg_set) - len(exp_set) - (bg_k - k)]]
    or_, p = fisher_exact(table, alternative='greater')
    results.append((term, k, bg_k, or_, p))
results.sort(key=lambda x: x[4])
if results:
    q = bh_fdr([r[4] for r in results])
    results = [r + (qq,) for r, qq in zip(results, q)]
# GO descriptions
go_desc = {}
for sp in ['pzi']:
    f = ANN / sp / 'swissprot_go.tsv'
    if f.exists():
        with f.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                go_desc.setdefault(r['go_id'], r['go_description'])
with (OUT / 'go_enrichment_expanded.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['go_id', 'go_description', 'expanded_families', 'bg_families',
                'odds_ratio', 'p', 'q'])
    for term, k, bg_k, or_, p, q in results:
        w.writerow([term, go_desc.get(term, ''), k, bg_k, f'{or_:.3f}', f'{p:.3g}', f'{q:.3g}'])
print('top GO terms:', [(r[0], go_desc.get(r[0], '')[:40], f'q={r[5]:.2g}') for r in results[:8]])
print('written to', OUT)
