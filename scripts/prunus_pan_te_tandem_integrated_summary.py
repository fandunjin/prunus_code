#!/usr/bin/env python3
"""Integrate standardized pan-TE composition with non-redundant TRF unions.

The output keeps two views separate but linked:

* broad RepeatMasker/pan-TE composition, where categories sum to the
  RepeatMasker union (this is the TE annotation view); and
* a mutually exclusive genome partition, where TRF and pan-TE overlap is
  assigned once to ``trf_pan_te_overlap`` (this is the genome accounting view).

This distinction is required because TRF-defined sequence and RepeatMasker
sequence overlap extensively in PZI.
"""

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
TE_CATEGORIES = ["LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA", "Helitron", "LINE", "SINE", "Other TE"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def num(row, field):
    try:
        return int(float(row.get(field, 0) or 0))
    except (TypeError, ValueError):
        return 0


def group_for(target):
    if target.startswith("pzi"): return "PZI"
    if target.startswith("psp"): return "PSP"
    if target.startswith("ppa"): return "PPA"
    if target == "pcam": return "PCAM"
    if target == "ppe": return "PPE"
    raise ValueError(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--te-composition", type=Path, required=True)
    parser.add_argument("--tandem-partition", type=Path, required=True)
    parser.add_argument("--tandem-period-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    te_rows = read_tsv(args.te_composition)
    partition_rows = read_tsv(args.tandem_partition)
    period_rows = read_tsv(args.tandem_period_summary)
    te_by_target = defaultdict(dict)
    for row in te_rows:
        te_by_target[row["target"]][row["category"]] = row
    partition = {row["target"]: row for row in partition_rows}
    missing = [target for target in TARGETS if target not in partition or target not in te_by_target]
    if missing:
        raise SystemExit("ERROR: missing integrated targets: " + ", ".join(missing))

    integrated = []
    for target in TARGETS:
        row = partition[target]
        genome_bp = num(row, "genome_bp")
        overlap_bp = num(row, "repeatmasker_overlap_bp")
        trf_bp = num(row, "trf_tandem_bp")
        te_only_bp = num(row, "pan_te_non_tandem_bp")
        neither_bp = num(row, "other_bp")
        trf_only_bp = trf_bp - overlap_bp
        te_total_bp = overlap_bp + te_only_bp
        if min(trf_only_bp, overlap_bp, te_only_bp, neither_bp) < 0:
            raise SystemExit(f"ERROR: negative integrated component for {target}")
        if overlap_bp + trf_only_bp + te_only_bp + neither_bp != genome_bp:
            raise SystemExit(f"ERROR: integrated partition does not sum to genome for {target}")
        out = {
            "target": target, "group": group_for(target), "genome_bp": genome_bp,
            "pan_te_total_bp": te_total_bp, "trf_total_bp": trf_bp,
            "trf_pan_te_overlap_bp": overlap_bp, "trf_only_bp": trf_only_bp,
            "pan_te_only_bp": te_only_bp, "neither_bp": neither_bp,
            "pan_te_total_fraction": f"{te_total_bp / genome_bp:.8f}",
            "trf_total_fraction": f"{trf_bp / genome_bp:.8f}",
            "trf_pan_te_overlap_fraction": f"{overlap_bp / genome_bp:.8f}",
            "trf_overlap_fraction": f"{overlap_bp / trf_bp:.8f}" if trf_bp else "0",
            "pan_te_overlap_fraction": f"{overlap_bp / te_total_bp:.8f}" if te_total_bp else "0",
        }
        for category in TE_CATEGORIES:
            out[f"te_{category.replace('/', '_').replace(' ', '_')}_bp"] = num(te_by_target[target].get(category, {}), "bp")
        integrated.append(out)

    group_rows = []
    for group in GROUPS:
        members = [row for row in integrated if row["group"] == group]
        out = {"group": group, "genome_packages": len(members)}
        for field in ("genome_bp", "pan_te_total_bp", "trf_total_bp", "trf_pan_te_overlap_bp", "trf_only_bp", "pan_te_only_bp", "neither_bp"):
            vals = [int(row[field]) for row in members]
            out[f"mean_{field}"] = round(statistics.mean(vals))
            out[f"min_{field}"] = min(vals)
            out[f"max_{field}"] = max(vals)
        for field in ("pan_te_total_fraction", "trf_total_fraction", "trf_pan_te_overlap_fraction", "trf_overlap_fraction", "pan_te_overlap_fraction"):
            vals = [float(row[field]) for row in members]
            out[f"mean_{field}"] = f"{statistics.mean(vals):.8f}"
            out[f"min_{field}"] = f"{min(vals):.8f}"
            out[f"max_{field}"] = f"{max(vals):.8f}"
        for category in TE_CATEGORIES:
            field = f"te_{category.replace('/', '_').replace(' ', '_')}_bp"
            out[f"mean_{field}"] = round(statistics.mean(int(row[field]) for row in members))
        group_rows.append(out)

    contrasts = []
    group_lookup = {row["group"]: row for row in group_rows}
    for group in ("PSP", "PPA", "PCAM", "PPE"):
        ref, comp = group_lookup["PZI"], group_lookup[group]
        fields = ("trf_pan_te_overlap_bp", "trf_only_bp", "pan_te_only_bp", "neither_bp")
        diffs = {field: int(ref[f"mean_{field}"]) - int(comp[f"mean_{field}"]) for field in fields}
        total = int(ref["mean_genome_bp"]) - int(comp["mean_genome_bp"])
        contrasts.append({
            "reference_group": "PZI", "comparator_group": group,
            "genome_size_difference_bp": total,
            "trf_pan_te_overlap_difference_bp": diffs["trf_pan_te_overlap_bp"],
            "trf_only_difference_bp": diffs["trf_only_bp"],
            "pan_te_only_difference_bp": diffs["pan_te_only_bp"],
            "neither_difference_bp": diffs["neither_bp"],
            "component_sum_bp": sum(diffs.values()),
            # Group means are rounded to integer bp before this contrast is
            # written.  Allow a small rounding residual rather than reporting
            # a false QC failure for a difference of only a few bases.
            "component_sum_matches_genome_difference": int(abs(sum(diffs.values()) - total) <= 10),
        })

    period_by_target = defaultdict(dict)
    for row in period_rows:
        period_by_target[row["target"]][row["period_class"]] = num(row, "exclusive_union_bp")
    period_out = []
    for group in GROUPS:
        members = [target for target in TARGETS if group_for(target) == group]
        periods = sorted(set().union(*(period_by_target[target] for target in members)), key=lambda x: ["1", "2-6", "7-20", "21-100", "101-500", "501-2000", ">2000"].index(x))
        for period in periods:
            period_out.append({"group": group, "period_class": period, "mean_exclusive_union_bp": round(statistics.mean(period_by_target[target].get(period, 0) for target in members))})

    fields = list(integrated[0])
    write_tsv(args.output_dir / "integrated_tandem_te_target.tsv", fields, integrated)
    write_tsv(args.output_dir / "integrated_tandem_te_group.tsv", list(group_rows[0]), group_rows)
    write_tsv(args.output_dir / "integrated_tandem_te_contrasts.tsv", list(contrasts[0]), contrasts)
    write_tsv(args.output_dir / "integrated_tandem_period_group.tsv", list(period_out[0]), period_out)
    manifest = {
        "te_composition_view": "RepeatMasker categories sum to the standardized RepeatMasker union and include TRF-overlapping sequence",
        "genome_accounting_view": "TRF-panTE overlap | TRF-only | panTE-only | neither; each base is assigned once",
        "warning": "pan-TE total and TRF total must not be added; the overlap is explicitly reported",
    }
    with open(args.output_dir / "integrated_tandem_te_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps({"targets": len(integrated), "groups": len(group_rows), "contrasts": len(contrasts), "output_dir": str(args.output_dir)}))


if __name__ == "__main__":
    main()
