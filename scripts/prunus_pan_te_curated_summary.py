#!/usr/bin/env python3
"""Apply the component-audit decisions to the pan-TE driver table.

The complete pan-library composition remains unchanged.  This command only
creates a conservative, audit-aware table for mechanistic TE-expansion claims.
"""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
PZI = TARGETS[:4]
PSP = TARGETS[4:8]
CONTEXT = TARGETS[8:]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def num(row, key):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return 0.0


def command(args):
    drivers = read_tsv(args.drivers)
    decisions = {r["pan_family"]: r for r in read_tsv(args.decisions) if r.get("pan_family")}
    missing = [r["pan_family"] for r in drivers if r["pan_family"] not in decisions]
    if missing and args.require_all_decisions:
        raise SystemExit("ERROR: no manual decision for: " + ", ".join(missing[:20]))

    intact = defaultdict(lambda: defaultdict(int))
    if args.intact.is_file():
        for r in read_tsv(args.intact):
            fam = r.get("pan_family", "")
            if fam:
                intact[fam][r.get("target", "")] += int(float(r.get("overlap_fraction", 0) or 0) >= args.intact_min_overlap)

    selected = []
    for r in drivers:
        fam = r["pan_family"]
        d = decisions.get(fam, {})
        decision = d.get("manual_decision", "")
        if decision != "include":
            continue
        out = dict(r)
        out["manual_decision"] = decision
        out["manual_evidence"] = d.get("manual_evidence", "")
        out["audit_notes"] = d.get("notes", "")
        out["pzi_intact_ltr_count"] = sum(intact[fam][t] for t in PZI)
        out["psp_intact_ltr_count"] = sum(intact[fam][t] for t in PSP)
        out["pzi_intact_ltr_consistent"] = int(all(intact[fam][t] > 0 for t in PZI))
        out["psp_intact_ltr_consistent"] = int(all(intact[fam][t] > 0 for t in PSP))
        selected.append(out)

    selected.sort(key=lambda r: num(r, "positive_excess_bp"), reverse=True)
    for i, r in enumerate(selected, 1):
        r["curated_rank"] = i
        r["curated_family"] = decisions[r["pan_family"]].get("curated_family", r["pan_family"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(drivers[0]) if drivers else []
    fields += ["manual_decision", "manual_evidence", "audit_notes", "curated_rank", "curated_family", "pzi_intact_ltr_count", "psp_intact_ltr_count", "pzi_intact_ltr_consistent", "psp_intact_ltr_consistent"]
    write_tsv(args.output_dir / "curated_te_driver_families.tsv", fields, selected)

    matrix_fields = ["pan_family", "curated_family", "class_family", "target", "exclusive_bp", "genome_fraction", "manual_evidence"]
    matrix = []
    for r in selected:
        for t in TARGETS:
            matrix.append({"pan_family": r["pan_family"], "curated_family": r["curated_family"], "class_family": r.get("class_family", ""), "target": t, "exclusive_bp": r.get(t + "_bp", "0"), "genome_fraction": "", "manual_evidence": r.get("manual_evidence", "")})
    write_tsv(args.output_dir / "curated_te_family_target_matrix.tsv", matrix_fields, matrix)

    by_class = defaultdict(lambda: {"families": 0, "positive_excess_bp": 0.0, "pzi_bp": 0.0, "psp_bp": 0.0, "context_bp": 0.0})
    for r in selected:
        c = r.get("class_family", "unknown")
        z = by_class[c]
        z["families"] += 1
        z["positive_excess_bp"] += num(r, "positive_excess_bp")
        z["pzi_bp"] += num(r, "pzi_mean_bp")
        z["psp_bp"] += num(r, "psp_mean_bp")
        z["context_bp"] += num(r, "ppa_pcam_ppe_median_bp")
    class_rows = [{"class_family": c, **{k: (int(v) if k == "families" else round(v)) for k, v in d.items()}} for c, d in sorted(by_class.items(), key=lambda x: x[1]["positive_excess_bp"], reverse=True)]
    write_tsv(args.output_dir / "curated_te_class_summary.tsv", ["class_family", "families", "positive_excess_bp", "pzi_bp", "psp_bp", "context_bp"], class_rows)

    manifest = {"drivers": str(args.drivers), "decisions": str(args.decisions), "intact_ltr": str(args.intact), "included_family_count": len(selected), "excluded_or_unreviewed_count": len(drivers) - len(selected), "policy": "only manual_decision=include is used for curated driver claims; full pan-library tables remain primary for composition"}
    with open(args.output_dir / "curated_te_manifest.json", "w", encoding="utf-8") as h:
        json.dump(manifest, h, indent=2, ensure_ascii=False)
    print(json.dumps(manifest, ensure_ascii=False))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--drivers", type=Path, required=True)
    p.add_argument("--decisions", type=Path, required=True)
    p.add_argument("--intact", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--intact-min-overlap", type=float, default=0.20)
    p.add_argument("--require-all-decisions", action="store_true")
    args = p.parse_args()
    for path in (args.drivers, args.decisions, args.intact):
        if not path.is_file():
            raise SystemExit(f"ERROR: missing {path}")
    command(args)


if __name__ == "__main__":
    main()
