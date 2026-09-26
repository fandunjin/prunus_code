#!/usr/bin/env python3
"""shared_specific_all17.py - clade set overlaps and their GO profiles with the
all-17-track annotation.

Sets (clade_family_scan_v3.tsv):
  R = racemose:  p_cs < 0.05, enr_cs > 1.3, n_species_up >= 2, te17 != True
  C = corymbose: p_cor_rest < 0.05, enr_cor_rest > 1.3
  S = solitary:  p_sol_rest < 0.05, enr_sol_rest > 1.3

Outputs (results_v2/):
  shared_specific_counts_all17.tsv   R_only/C_only/S_only/R_C/R_S/C_S/R_C_S
  shared_specific_go_all17.tsv       per-group GO enrichment (Fisher + BH)
"""
import csv
from collections import Counter
from pathlib import Path

from scipy.stats import fisher_exact

_CANDIDATES = [Path(r'${PROJECT_ROOT}'), Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'
OUT = B / 'analysis/racemose_inflorescence_scan_20260922/results_v2'
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


def main():
    rows = list(csv.DictReader((FIG / 'data/clade_family_scan_v3.tsv').open(), delimiter='\t'))
    fam_go = {r['og']: set(x for x in (r['go17'] or '').split(';') if x) for r in rows}

    R = {r['og'] for r in rows if num(r, 'q_cs') < 0.05 and num(r, 'enr_cs') > 1.3
         and int(r['n_species_up']) >= 2 and r['te17'] != 'True'}
    C = {r['og'] for r in rows if num(r, 'q_cor_rest') < 0.05 and num(r, 'enr_cor_rest') > 1.3}
    S = {r['og'] for r in rows if num(r, 'q_sol_rest') < 0.05 and num(r, 'enr_sol_rest') > 1.3}

    groups = {
        'R_only': R - C - S,
        'C_only': C - R - S,
        'S_only': S - R - C,
        'C_S': (C & S) - R,
    }
    counts = {
        'R_only': len(groups['R_only']), 'C_only': len(groups['C_only']),
        'S_only': len(groups['S_only']),
        'R_C': len((R & C) - S), 'R_S': len((R & S) - C),
        'C_S': len((C & S) - R), 'R_C_S': len(R & C & S),
    }
    print('counts:', counts, '| totals R/C/S:', len(R), len(C), len(S))

    with (OUT / 'shared_specific_counts_all17.tsv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['category', 'n_families'])
        for k, v in counts.items():
            w.writerow([k, v])

    go_desc = {}
    if GO_DESC_SRC.exists():
        with GO_DESC_SRC.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                go_desc.setdefault(r['go_id'], r['go_description'])

    bg_set = {og for og, gos in fam_go.items() if gos}
    out_rows = []
    for label in ('R_only', 'C_only', 'S_only', 'C_S'):
        members = groups[label] & bg_set
        counts_go = Counter()
        for og in members:
            for g in fam_go[og]:
                counts_go[g] += 1
        results = []
        for term, k in counts_go.items():
            if k < 3:
                continue
            bg_k = sum(1 for og in bg_set if term in fam_go[og])
            table = [[k, len(members) - k], [bg_k - k, len(bg_set) - len(members) - (bg_k - k)]]
            _, p = fisher_exact(table, alternative='greater')
            results.append((term, k, bg_k, p))
        results.sort(key=lambda x: x[3])
        q = bh_fdr([r[3] for r in results]) if results else []
        for (term, k, bg_k, p), qq in zip(results, q):
            out_rows.append([label, term, go_desc.get(term, ''), k, bg_k, p, qq])
        top = [(t, go_desc.get(t, '')[:38], k, bg_k, f'{qq:.1e}')
               for (t, k, bg_k, _), qq in zip(results[:6], q[:6])]
        print(f'{label}: n={len(members)}, terms={len(results)}')
        for t in top:
            print('   ', t)

    with (OUT / 'shared_specific_go_all17.tsv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['group', 'go_id', 'go_description', 'k', 'bg_k', 'p', 'q'])
        for r in out_rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4], f'{r[5]:.3g}', f'{r[6]:.3g}'])
    print('written shared_specific_counts_all17.tsv / shared_specific_go_all17.tsv')


if __name__ == '__main__':
    main()
