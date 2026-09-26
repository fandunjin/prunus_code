#!/usr/bin/env python3
"""Aggregate independently estimated CodeML Ks values for five Figure 2 panels."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


PANEL_MAP = {
    "ab_paralogue": "ab_paralogue",
    "orthologue_pcam": "orthologue_pcam",
    "orthologue_ppe": "orthologue_ppe",
    "orthologue_gtr": "orthologue_gtr",
    "self_paralogue": "self_paralogue",
}


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("Missing header: {}".format(path))
        return list(reader)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError("Missing input: {}".format(path))


def manifests(input_root):
    rows = []
    for group in ("ab", "external", "self"):
        manifest = input_root / group / "comparisons.tsv"
        require_file(manifest)
        values = read_tsv(manifest)
        for row in values:
            if row.get("group") != group:
                raise RuntimeError("Group mismatch in {}".format(manifest))
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    output_rows = []
    audit_rows = []
    seen = set()
    for planned in manifests(args.input_root):
        comparison = planned["comparison"]
        group = planned["group"]
        result = args.results_root / group / comparison
        complete = result / "COMPLETE"
        retained = result / "audit" / (comparison + ".ks_0_5.tsv")
        parameters = result / "run_parameters.tsv"
        require_file(complete)
        require_file(retained)
        require_file(parameters)
        parameter_text = parameters.read_text(encoding="utf-8", errors="replace")
        for expected in ("MUSCLE_3.8.1551", "codeml_CodonFreq\t2_F3X4", "pairwise\ttrue"):
            if expected not in parameter_text:
                raise RuntimeError("CodeML contract missing from {}: {}".format(parameters, expected))
        source_rows = read_tsv(retained)
        panel = PANEL_MAP.get(planned.get("panel", ""))
        curve = planned.get("curve", "")
        if not panel or not curve:
            raise RuntimeError("Invalid plotting metadata for {}".format(comparison))
        count = 0
        for row in source_rows:
            try:
                ks = float(row.get("dS", ""))
            except ValueError:
                raise RuntimeError("Invalid dS in {}".format(retained))
            if not math.isfinite(ks) or not (0.0 < ks < 5.0):
                raise RuntimeError("Unfiltered Ks in {}: {}".format(retained, ks))
            family_id = row.get("family_id", "")
            if not family_id:
                raise RuntimeError("Missing family_id in {}".format(retained))
            unique = (comparison, family_id)
            if unique in seen:
                raise RuntimeError("Duplicate combined pair family: {}".format(unique))
            seen.add(unique)
            output_rows.append({
                "panel": panel,
                "curve": curve,
                "source_comparison": comparison,
                "family_id": family_id,
                "query_label": planned.get("query_label", ""),
                "ref_label": planned.get("ref_label", ""),
                "gene_1": row.get("gene_1", ""),
                "gene_2": row.get("gene_2", ""),
                "Ks": "{:.10g}".format(ks),
            })
            count += 1
        if count == 0:
            raise RuntimeError("No retained CodeML Ks rows for {}".format(comparison))
        audit_rows.append({
            "group": group,
            "panel": panel,
            "curve": curve,
            "source_comparison": comparison,
            "retained_0_lt_Ks_lt_5": count,
            "source_table": str(retained),
        })

    if not output_rows:
        raise SystemExit("[ERROR] No CodeML Ks values were collected")
    output_rows.sort(key=lambda row: (row["panel"], row["curve"], row["source_comparison"], row["family_id"]))
    write_tsv(
        args.outdir / "five_panel_ks.tsv",
        ["panel", "curve", "source_comparison", "family_id", "query_label", "ref_label", "gene_1", "gene_2", "Ks"],
        output_rows,
    )
    write_tsv(
        args.outdir / "comparison_audit.tsv",
        ["group", "panel", "curve", "source_comparison", "retained_0_lt_Ks_lt_5", "source_table"],
        audit_rows,
    )
    curve_counts = defaultdict(int)
    for row in output_rows:
        curve_counts[(row["panel"], row["curve"])] += 1
    write_tsv(
        args.outdir / "curve_pair_counts.tsv",
        ["panel", "curve", "retained_0_lt_Ks_lt_5"],
        [{"panel": key[0], "curve": key[1], "retained_0_lt_Ks_lt_5": value}
         for key, value in sorted(curve_counts.items())],
    )
    (args.outdir / "method_contract.txt").write_text(
        "All 48 phased comparisons were estimated separately with the old Figure 2c CodeML contract.\n"
        "Alignment: MUSCLE 3.8.1551 protein alignment followed by PAL2NAL codon backtranslation.\n"
        "Ks estimator: PAML CodeML pairwise runmode=-2, CodonFreq=2 (F3X4), cleandata=1.\n"
        "Pair sources: corrected JCVI anchors for external/A-B pairs and i-ADHoRe cross-chromosome anchorpoints for self pairs.\n"
        "Filter: 0<Ks<5. PZI and PSP A1/A2 or B1/B2 values are pooled only after independent pair-level estimation.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("CodeML Figure 2 five-panel data combined.\n", encoding="ascii")
    print("[OK] Combined {} CodeML Ks pairs from {} planned comparisons: {}".format(
        len(output_rows), len(audit_rows), args.outdir))


if __name__ == "__main__":
    main()
