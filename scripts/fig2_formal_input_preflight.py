#!/usr/bin/env python3
"""Validate the data contract for the formal Figure 2a-d renderer.

This checker intentionally treats anchor-state summaries and candidate joint
states as insufficient for panels c and d.  Formal c requires FraBias window
retention output; formal d requires POInT likelihood/pillar output.
"""

import argparse
import csv
from pathlib import Path


def read_header(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"), [])


def check(path, required):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return "MISSING", []
    header = read_header(p)
    missing = [name for name in required if name not in header]
    return ("PASS" if not missing else "FAIL_COLUMNS"), missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", required=True)
    ap.add_argument("--topologies", required=True)
    ap.add_argument("--support", required=True)
    ap.add_argument("--pug", required=True)
    ap.add_argument("--frabias", required=True)
    ap.add_argument("--point", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    checks = [
        ("a_ks", args.ks, ["comparison", "Ks", "kaks_status"]),
        ("b_topologies", args.topologies, ["contrast", "topology", "frequency"]),
        ("b_support", args.support, ["displayed_descendant_taxa", "astral_pp", "gcf", "likelihood_scf"]),
        ("b_pug", args.pug, ["event", "node", "support_cutoff", "gene_tree_count"]),
        ("c_frabias", args.frabias, ["pair_label", "track_id", "target_seqid", "window_index", "retention_proportion"]),
        ("d_point", args.point, ["branch_id", "subgenome", "surviving_genes", "losses", "retention_ratio", "posterior", "analysis_status"]),
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    formal = True
    lines = ["Figure 2 formal input preflight", ""]
    for name, path, required in checks:
        status, missing = check(path, required)
        lowered = str(path).lower()
        if name == "a_ks" and any(token in lowered for token in ("jcvi", "kaks_calculator", "all_pairs_yn", "collinear_all_pairs_yn")):
            status, missing = "FAIL_FORBIDDEN_SOURCE", ["wgd/CodeML/i-ADHoRe Ks output"]
        if name == "b_pug" and "iadhore" not in lowered:
            status, missing = "FAIL_FORBIDDEN_SOURCE", ["PUG mapping from i-ADHoRe collinear pairs"]
        if name == "c_frabias" and any(token in lowered for token in ("jcvi", "candidate", "alltaxa_fractionation_states")):
            status, missing = "FAIL_FORBIDDEN_SOURCE", ["formal LAST/Quota/FraBias path"]
        if name == "d_point" and any(token in lowered for token in ("jcvi", "dollo", "candidate", "alltaxa_fractionation_states")):
            status, missing = "FAIL_FORBIDDEN_SOURCE", ["official converged POInT path"]
        if status != "PASS":
            formal = False
        lines.append("{}\t{}\t{}".format(name, status, ",".join(missing) or path))
    lines += [
        "",
        "formal_figure_ready\t{}".format("YES" if formal else "NO"),
        "panel_c_boundary\tanchor-state summaries are not FraBias retention output",
        "panel_d_boundary\tcandidate joint states are not POInT likelihood output",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("[{}] Formal Figure 2 input preflight: {}".format("OK" if formal else "BLOCKED", out))
    if not formal:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
