#!/usr/bin/env python3
"""dnds_branch_permutation_control.py - permuted-branch calibration control for the
racemose branch dN/dS test (R8, 2026-09-26).

Design: the same random single-copy orthogroups and codon alignments as the
genome-wide background control are re-tested with the foreground labels (#1)
moved from the racemose clade to every subset of size 3-5 of the five
non-racemose branches. The null distribution of the fraction of sets with
p < 0.05 and of the foreground/background omega ratio is compared with the
observed racemose values.

Usage (WSL):
  python3 dnds_branch_permutation_control.py --bg B --workers 8
"""
import argparse
import csv
import itertools
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dnds_background_controls as bg  # noqa: E402  (main() is guarded)

OUT = Path('${PROJECT_ROOT}/analysis/dnds_branch_20260922')

RACEMOSE = {'ppa', 'psp', 'pzi', 'pspzi', 'ppa_pspzi'}
ALTERNATIVES = ['gtr', 'pcam', 'ppe', 'pcamppe', 'root_stem']


def tree_with_foreground(slots):
    def L(name, text):
        return text + ' #1' if name in slots else text

    psp = L('psp', 'psp:0.1')
    pzi = L('pzi', 'pzi:0.1')
    pspzi = L('pspzi', '({}, {}):0.05'.format(psp, pzi))
    ppa = L('ppa', 'ppa:0.1')
    rac = L('ppa_pspzi', '({}, {}):0.05'.format(ppa, pspzi))
    pcam = L('pcam', 'pcam:0.1')
    ppe = L('ppe', 'ppe:0.1')
    pcamppe = L('pcamppe', '({}, {}):0.1'.format(pcam, ppe))
    root = L('root_stem', '({}, {}):0.2'.format(rac, pcamppe))
    gtr = L('gtr', 'gtr:0.3')
    return '({}, {});'.format(gtr, root)


def run_one(job):
    run_dir, aln_phy, tree_text, model, tag = job
    return bg.run_codeml(Path(run_dir), Path(aln_phy), tree_text, model, tag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bg', default='B', choices=['A', 'B'])
    ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()

    runs_root = OUT / ('runs_bg' if args.bg == 'A' else 'runs_bg_B')
    og_dirs = sorted([d for d in runs_root.iterdir() if (d / 'aln.phy').is_file()])
    print('background {}: {} alignments'.format(args.bg, len(og_dirs)))

    jobs = [(str(d), str(d / 'aln.phy'), bg.TREE6_NULL, 0, 'null') for d in og_dirs]
    null_lnl = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for d, res in zip(og_dirs, ex.map(run_one, jobs)):
            if res:
                null_lnl[d.name] = res['lnL']
    print('null fits:', len(null_lnl))

    subsets = []
    for k in (3, 4, 5):
        subsets.extend(itertools.combinations(ALTERNATIVES, k))
    print('permuted foregrounds:', len(subsets))

    results = []
    for i, slots in enumerate(subsets):
        tree_text = tree_with_foreground(set(slots))
        jobs = [(str(d), str(d / 'aln.phy'), tree_text, 2, 'perm{}'.format(i))
                for d in og_dirs if d.name in null_lnl]
        names = [d for d in og_dirs if d.name in null_lnl]
        ps, ratios = [], []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for d, res in zip(names, ex.map(run_one, jobs)):
                if not res:
                    continue
                dlnl = 2 * (res['lnL'] - null_lnl[d.name])
                ps.append(float(chi2.sf(max(dlnl, 0), 1)))
                if res['omega_bg']:
                    ratios.append(res['omega_fg'] / res['omega_bg'])
        if ps:
            results.append(dict(perm=i, slots='+'.join(slots), n=len(ps),
                                frac_sig=sum(1 for p in ps if p < 0.05) / len(ps),
                                mean_ratio=sum(ratios) / len(ratios) if ratios else float('nan')))
            print('perm {} [{}]: n={} sig={:.1%} mean_ratio={:.2f}'.format(
                i, results[-1]['slots'], len(ps), results[-1]['frac_sig'], results[-1]['mean_ratio']), flush=True)

    out = OUT / 'results/branch_dnds_permutation_control_{}.tsv'.format(args.bg)
    with out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    fs = sorted(r['frac_sig'] for r in results)
    mr = sorted(r['mean_ratio'] for r in results)
    print('null distribution ({} permutations):'.format(len(results)))
    print('  fraction significant: min {:.1%} median {:.1%} max {:.1%}'.format(fs[0], fs[len(fs) // 2], fs[-1]))
    print('  mean omega ratio:     min {:.2f} median {:.2f} max {:.2f}'.format(mr[0], mr[len(mr) // 2], mr[-1]))
    print('written', out)


if __name__ == '__main__':
    main()
