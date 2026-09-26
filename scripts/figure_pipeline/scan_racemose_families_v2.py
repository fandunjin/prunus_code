#!/usr/bin/env python3
"""scan_racemose_families_v2.py - refined genome-wide clade scan.

Improvements over v1
  * per-track density null (each haplotype track ~ equal gene number), so the
    uniform WGD doubling is the baseline and deviations flag family-specific
    retention/fractionation;
  * per-species consistency: the racemose test must be significant in >=2 of
    PZI/PSP/PPA (not just in the pooled sum), removing single-species artefacts;
  * TE-derived families are flagged (transposon keywords) and reported separately;
  * symmetric tests for the corymbose and solitary clades against the rest.

Groups
  racemose : ppa_A, ppa_B, psp_A, psp_B, pzi_A, pzi_B   (A/B halves also tracked)
  corymbose: pcam, pav, pspe, pmal
  solitary : par, pmu, pdu, psal, ppe
  outgroup : gtr, rch

Outputs (results_v2/)
  clade_family_scan.tsv
  racemose_expanded_consistent.tsv
  racemose_te_families.tsv
  clade_specific_expansions.tsv
  go_enrichment_consistent.tsv
"""
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import binomtest, fisher_exact

B = Path('${PROJECT_ROOT}')
ORTHO = B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups'
FROZEN = B / '0903/analysisphylogenetics/01_frozen_inputs/main/id_maps'
ANN = B / '0903/annotation/functional'
OUT = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'
OUT.mkdir(parents=True, exist_ok=True)

TR = {
    'RAC': ['ppa_A', 'ppa_B', 'psp_A', 'psp_B', 'pzi_A', 'pzi_B'],
    'RAC_A': ['ppa_A', 'psp_A', 'pzi_A'],
    'RAC_B': ['ppa_B', 'psp_B', 'pzi_B'],
    'PZI': ['pzi_A', 'pzi_B'],
    'PSP': ['psp_A', 'psp_B'],
    'PPA': ['ppa_A', 'ppa_B'],
    'COR': ['pcam', 'pav', 'pspe', 'pmal'],
    'SOL': ['par', 'pmu', 'pdu', 'psal', 'ppe'],
    'OUT': ['gtr', 'rch'],
}
TE_RE = re.compile(r'transposon|retrotransposon|retrovirus|reverse transcriptase|'
                   r'integrase|transposase|gag-polyprotein|copia|gypsy|'
                   r'pol polyprotein|rnase h', re.I)


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        prev = min(prev, pvals[idx] * n / i)
        q[idx] = prev
    return q


def load_id_map(track):
    m = {}
    f = FROZEN / f'{track}.id_map.tsv'
    if f.exists():
        with f.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                m[r['stable_id']] = r['original_id']
    return m


IDMAP = {t: load_id_map(t) for t in TR['RAC']}
DESC = {}
for sp in ['pzi', 'psp', 'ppa']:
    d = {}
    f = ANN / sp / 'swissprot_functional_annotation.tsv'
    if f.exists():
        with f.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                d[r['protein_id']] = {
                    'desc': (r.get('protein_description') or '').split(';')[0][:150],
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
    for tr in TR['RAC']:
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
            ' | '.join(dict.fromkeys(descs))[:150], gos)


# ---------------------------------------------------------------- counts
counts, totals = {}, Counter()
with (ORTHO / 'Orthogroups.GeneCount.tsv').open(errors='ignore') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        og = r['Orthogroup']
        d = {t: int(r[t]) for t in sum(TR.values(), [])}
        # keep only unique tracks
        d = {k: v for k, v in d.items()}
        counts[og] = d
        for t, v in d.items():
            totals[t] += v

TOT = {g: sum(totals[t] for t in tracks) for g, tracks in TR.items()}
print('group gene totals:', {g: TOT[g] for g in ('RAC', 'COR', 'SOL', 'OUT')})


def test(obs, n, base, base_tot):
    if base_tot == 0 or n == 0:
        return None
    p = min(1.0, base / base_tot)
    exp = p * n
    return exp, binomtest(obs, n, p, alternative='greater').pvalue, \
        binomtest(obs, n, p, alternative='less').pvalue


rows = []
for og, d in counts.items():
    s = {g: sum(d[t] for t in tracks) for g, tracks in TR.items()}
    present = sum(1 for t in sum(TR.values(), []) if d[t] > 0)
    if present < 2 or s['RAC'] + s['COR'] + s['SOL'] + s['OUT'] > 150:
        continue
    if s['COR'] + s['SOL'] + s['OUT'] < 3:
        continue
    r = {}
    for g, base_groups in [('vs_cor_sol', ('COR', 'SOL')), ('vs_cor', ('COR',)),
                           ('vs_sol', ('SOL',)), ('vs_out', ('OUT',))]:
        base = sum(s[x] for x in base_groups)
        base_tot = sum(TOT[x] for x in base_groups)
        r[g] = test(s['RAC'], TOT['RAC'], base, base_tot)
    # per-species consistency
    cons = {}
    for sp in ('PZI', 'PSP', 'PPA'):
        base = s['COR'] + s['SOL']
        cons[sp] = test(s[sp], TOT[sp], base, TOT['COR'] + TOT['SOL'])
    # symmetric clade tests (corymbose / solitary against the rest per track)
    rest_for_cor = ('RAC', 'SOL', 'OUT')
    r['cor_vs_rest'] = test(s['COR'], TOT['COR'],
                            sum(s[x] for x in rest_for_cor),
                            sum(TOT[x] for x in rest_for_cor))
    rest_for_sol = ('RAC', 'COR', 'OUT')
    r['sol_vs_rest'] = test(s['SOL'], TOT['SOL'],
                            sum(s[x] for x in rest_for_sol),
                            sum(TOT[x] for x in rest_for_sol))
    rows.append((og, s, r, cons))
print('families tested:', len(rows))

# annotate screened families
screen = [x for x in rows if
          (x[2]['vs_cor_sol'] and x[2]['vs_cor_sol'][1] < 0.05 and
           x[2]['vs_cor_sol'][0] / (x[1]['RAC'] or 1) > 0) or
          (x[2]['vs_cor_sol'] and x[2]['vs_cor_sol'][2] < 0.05) or
          (x[2]['cor_vs_rest'] and x[2]['cor_vs_rest'][1] < 0.05) or
          (x[2]['sol_vs_rest'] and x[2]['sol_vs_rest'][1] < 0.05)]
print('screened for annotation:', len(screen))
annot = {og: annotate_og(og) for og, *_ in screen}


def enr(t):
    if not t or t[0] in (0, float('nan')):
        return float('nan')
    return t[0] / t[0]  # placeholder, replaced below


def enrich(obs, n, base, base_tot):
    if base_tot == 0 or n == 0 or base == 0:
        return float('nan')
    return (obs / n) / (base / base_tot)


def pack(og, s, r, cons, sym, desc, te):
    e_cs = enrich(s['RAC'], TOT['RAC'], s['COR'] + s['SOL'], TOT['COR'] + TOT['SOL'])
    e_cor = enrich(s['RAC'], TOT['RAC'], s['COR'], TOT['COR'])
    e_sol = enrich(s['RAC'], TOT['RAC'], s['SOL'], TOT['SOL'])
    e_out = enrich(s['RAC'], TOT['RAC'], s['OUT'], TOT['OUT'])
    n_sig = sum(1 for sp in ('PZI', 'PSP', 'PPA')
                if cons[sp] and cons[sp][1] < 0.05 and enrich(s[sp], TOT[sp], s['COR'] + s['SOL'],
                                                              TOT['COR'] + TOT['SOL']) > 1)
    return dict(og=og, rac=s['RAC'], racA=s['RAC_A'], racB=s['RAC_B'],
                pzi=s['PZI'], psp=s['PSP'], ppa=s['PPA'],
                cor=s['COR'], sol=s['SOL'], out=s['OUT'],
                enr_cs=e_cs, p_cs=r['vs_cor_sol'][1] if r['vs_cor_sol'] else 1,
                p_cs_less=r['vs_cor_sol'][2] if r['vs_cor_sol'] else 1,
                enr_cor=e_cor, p_cor=r['vs_cor'][1] if r['vs_cor'] else 1,
                enr_sol=e_sol, p_sol=r['vs_sol'][1] if r['vs_sol'] else 1,
                enr_out=e_out, p_out=r['vs_out'][1] if r['vs_out'] else 1,
                n_species_up=n_sig,
                enr_cor_rest=(enrich(s['COR'], TOT['COR'],
                                     sum(s[x] for x in ('RAC', 'SOL', 'OUT')),
                                     sum(TOT[x] for x in ('RAC', 'SOL', 'OUT')))),
                p_cor_rest=r['cor_vs_rest'][1] if r['cor_vs_rest'] else 1,
                enr_sol_rest=(enrich(s['SOL'], TOT['SOL'],
                                     sum(s[x] for x in ('RAC', 'COR', 'OUT')),
                                     sum(TOT[x] for x in ('RAC', 'COR', 'OUT')))),
                p_sol_rest=r['sol_vs_rest'][1] if r['sol_vs_rest'] else 1,
                symbol=sym, description=desc, te=te)


recs = []
for og, s, r, cons in rows:
    sym, desc, gos = annot.get(og, ('', '', set()))
    te = bool(TE_RE.search(sym + ' ' + desc)) if (sym or desc) else False
    d = pack(og, s, r, cons, sym, desc, te)
    recs.append((d, gos))

# Benjamini-Hochberg FDR across all tested families (2026-09-26: multiplicity control)
q_cs = bh_fdr([d['p_cs'] for d, _ in recs])
q_cs_less = bh_fdr([d['p_cs_less'] for d, _ in recs])
q_cor_rest = bh_fdr([d['p_cor_rest'] for d, _ in recs])
q_sol_rest = bh_fdr([d['p_sol_rest'] for d, _ in recs])
for (d, _), qv, qv2, qc, qs in zip(recs, q_cs, q_cs_less, q_cor_rest, q_sol_rest):
    d['q_cs'] = qv
    d['q_cs_less'] = qv2
    d['q_cor_rest'] = qc
    d['q_sol_rest'] = qs

cols = ['og', 'rac', 'racA', 'racB', 'pzi', 'psp', 'ppa', 'cor', 'sol', 'out',
        'enr_cs', 'p_cs', 'q_cs', 'p_cs_less', 'q_cs_less', 'enr_cor', 'p_cor', 'enr_sol', 'p_sol', 'enr_out', 'p_out',
        'n_species_up', 'enr_cor_rest', 'p_cor_rest', 'q_cor_rest', 'enr_sol_rest', 'p_sol_rest', 'q_sol_rest',
        'symbol', 'description', 'te']
with (OUT / 'clade_family_scan.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(cols)
    for d, _ in recs:
        w.writerow([f'{d[c]:.4g}' if isinstance(d[c], float) else d[c] for c in cols])

exp_consistent = [d for d, _ in recs
                  if d['q_cs'] < 0.05 and d['enr_cs'] > 1.3 and d['n_species_up'] >= 2 and not d['te']]
exp_te = [d for d, _ in recs if d['te'] and d['q_cs'] < 0.05 and d['enr_cs'] > 1.3]
cor_spec = [d for d, _ in recs if d['q_cor_rest'] < 0.05 and d['enr_cor_rest'] > 1.3]
sol_spec = [d for d, _ in recs if d['q_sol_rest'] < 0.05 and d['enr_sol_rest'] > 1.3]
print(f'consistent racemose-expanded (non-TE): {len(exp_consistent)}; TE families: {len(exp_te)}; '
      f'corymbose-specific: {len(cor_spec)}; solitary-specific: {len(sol_spec)}')


def dump(name, items):
    with (OUT / name).open('w', newline='') as fh:
        cols = ['og', 'rac', 'racA', 'racB', 'pzi', 'psp', 'ppa', 'cor', 'sol', 'out',
                'enr_cs', 'p_cs', 'q_cs', 'p_cs_less', 'q_cs_less', 'enr_cor', 'p_cor', 'enr_sol', 'p_sol', 'enr_out',
                'n_species_up', 'enr_cor_rest', 'p_cor_rest', 'q_cor_rest',
                'enr_sol_rest', 'p_sol_rest', 'q_sol_rest', 'symbol', 'description']
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for d in sorted(items, key=lambda x: -x['enr_cs']):
            w.writerow(d)


dump('racemose_expanded_consistent.tsv', exp_consistent)
dump('racemose_te_families.tsv', exp_te)
dump('corymbose_specific.tsv', cor_spec)
dump('solitary_specific.tsv', sol_spec)

# ---------------------------------------------------------------- GO enrichment
fam_go = {d['og']: gos for d, gos in recs if gos}
exp_set = {d['og'] for d in exp_consistent if d['og'] in fam_go}
bg_set = set(fam_go)
print(f'GO families: {len(bg_set)}; expanded with GO: {len(exp_set)}')
go_count = Counter()
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
go_desc = {}
f = ANN / 'pzi/swissprot_go.tsv'
if f.exists():
    with f.open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            go_desc.setdefault(r['go_id'], r['go_description'])
with (OUT / 'go_enrichment_consistent.tsv').open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['go_id', 'go_description', 'expanded', 'bg', 'odds_ratio', 'p', 'q'])
    for term, k, bg_k, or_, p, qq in results:
        w.writerow([term, go_desc.get(term, ''), k, bg_k, f'{or_:.3f}', f'{p:.3g}', f'{qq:.3g}'])
print('top GO:', [(r[0], go_desc.get(r[0], '')[:35], f'q={r[5]:.2g}') for r in results[:10]])
print('written', OUT)
