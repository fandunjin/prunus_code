#!/usr/bin/env python3
"""go_enrichment_clades_all17.py - GO enrichment for corymbose / solitary / shared
sets using the all-17-track annotation.

Sets (from clade_family_scan_v3.tsv):
  corymbose-specific : p_cor_rest < 0.05 and enr_cor_rest > 1.3
  solitary-specific  : p_sol_rest < 0.05 and enr_sol_rest > 1.3
  shared (cor & sol) : intersection of the two
  racemose-consistent: p_cs < 0.05, enr_cs > 1.3, n_species_up >= 2 (non-TE) - recomputed
                       with the all-17 annotation for comparison

Outputs (analysis/racemose_inflorescence_scan_20260922/results_v2/):
  go_enrichment_corymbose_all17.tsv
  go_enrichment_solitary_all17.tsv
  go_enrichment_cor_sol_shared_all17.tsv
  go_enrichment_racemose_all17.tsv
  cor_sol_shared_annotated.tsv
"""
import csv
import math
from collections import Counter
from pathlib import Path

from scipy.stats import fisher_exact

_CANDIDATES = [Path(r'${PROJECT_ROOT}'), Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'
SCAN = B / 'analysis/racemose_inflorescence_scan_20260922'
OUT = SCAN / 'results_v2'
GO_DESC_SRC = B / '0903/annotation/functional/pzi/swissprot_go.tsv'


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


def num(r, k):
    try:
        return float(r[k])
    except Exception:
        return float('nan')


def go_sets(rows):
    return {r['og']: set(x for x in (r['go17'] or '').split(';') if x) for r in rows}


def enrich(rows, set_ogs, fam_go, go_desc, label):
    bg_ogs = [og for og in fam_go if fam_go[og]]
    bg_set = set(bg_ogs)
    set_ogs = set(set_ogs) & bg_set
    counts = Counter()
    for og in set_ogs:
        for g in fam_go[og]:
            counts[g] += 1
    results = []
    for term, k in counts.items():
        if k < 3:
            continue
        bg_k = sum(1 for og in bg_set if term in fam_go[og])
        table = [[k, len(set_ogs) - k], [bg_k - k, len(bg_set) - len(set_ogs) - (bg_k - k)]]
        _, p = fisher_exact(table, alternative='greater')
        results.append((term, k, bg_k, p))
    results.sort(key=lambda x: x[3])
    if results:
        q = bh_fdr([r[3] for r in results])
    else:
        q = []
    with (OUT / f'go_enrichment_{label}_all17.tsv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['go_id', 'go_description', 'expanded', 'bg', 'p', 'q'])
        for (term, k, bg_k, p), qq in zip(results, q):
            w.writerow([term, go_desc.get(term, ''), k, bg_k, f'{p:.3g}', f'{qq:.3g}'])
    top = [(r[0], go_desc.get(r[0], '')[:42], r[1], r[2], f'{qq:.1e}')
           for r, qq in zip(results[:10], q[:10])]
    print(f'{label}: set={len(set_ogs)}, background={len(bg_set)}, terms(k>=3)={len(results)}')
    for t in top:
        print('   ', t)
    return results, q


def main():
    rows = list(csv.DictReader((FIG / 'data/clade_family_scan_v3.tsv').open(), delimiter='\t'))
    fam_go = go_sets(rows)

    go_desc = {}
    if GO_DESC_SRC.exists():
        with GO_DESC_SRC.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                go_desc.setdefault(r['go_id'], r['go_description'])

    cor_spec = [r['og'] for r in rows
                if num(r, 'p_cor_rest') < 0.05 and num(r, 'enr_cor_rest') > 1.3]
    sol_spec = [r['og'] for r in rows
                if num(r, 'p_sol_rest') < 0.05 and num(r, 'enr_sol_rest') > 1.3]
    shared = sorted(set(cor_spec) & set(sol_spec))
    rac_set = [r['og'] for r in rows
               if num(r, 'p_cs') < 0.05 and num(r, 'enr_cs') > 1.3
               and int(r['n_species_up']) >= 2 and r['te17'] != 'True']

    print(f'corymbose {len(cor_spec)}, solitary {len(sol_spec)}, shared {len(shared)}, '
          f'racemose {len(rac_set)}')
    enrich(rows, cor_spec, fam_go, go_desc, 'corymbose')
    enrich(rows, sol_spec, fam_go, go_desc, 'solitary')
    enrich(rows, shared, fam_go, go_desc, 'cor_sol_shared')
    enrich(rows, rac_set, fam_go, go_desc, 'racemose')

    # annotated shared table
    by_og = {r['og']: r for r in rows}
    shared_rows = [by_og[og] for og in shared]
    shared_rows.sort(key=lambda r: -min(num(r, 'enr_cor_rest'), num(r, 'enr_sol_rest')))
    with (OUT / 'cor_sol_shared_annotated.tsv').open('w', newline='', encoding='utf-8') as fh:
        cols = ['og', 'symbol17', 'description17', 'te17', 'pzi', 'psp', 'ppa',
                'cor', 'sol', 'out', 'enr_cor_rest', 'enr_sol_rest', 'p_cor_rest',
                'p_sol_rest', 'go17', 'src_tracks17']
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in shared_rows:
            w.writerow(r)
    named = sum(1 for r in shared_rows if (r['symbol17'] or '').strip())
    print(f'cor_sol_shared_annotated.tsv: {len(shared_rows)} families, '
          f'{named} with symbol, '
          f'{sum(1 for r in shared_rows if r["go17"])} with GO')


if __name__ == '__main__':
    main()
