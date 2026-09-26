#!/usr/bin/env python3
"""Prepare auditable inputs for the three-panel, piecewise-axis Ks figure.

Panel a keeps only the existing representative A-B homoeologous pairs.
Panel b reuses the completed three-reference YN Ks table without changing a
single pair. Panel c reads six independently run subgenome self-synteny Ks
analyses, rather than subsetting the earlier composite A+B analyses.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


AB_LABELS = ("PZI", "PSP", "PPA")
SUBGENOME_SELF = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B")
SINGLE_SELF = ("PCAM", "PPE", "GTR")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Missing header in {}".format(path))
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_ks(value: str, path: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise ValueError("Invalid Ks in {}: {}".format(path, value))
    if not math.isfinite(result) or not 0.0 < result < 5.0:
        raise ValueError("Expected 0<Ks<5 in {} but found {}".format(path, value))
    return result


def panel_a_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(path):
        label = row.get("comparison", "")
        if label not in AB_LABELS:
            continue
        value = parse_ks(row.get("Ks", ""), path)
        rows.append({
            "comparison": label,
            "family_id": row.get("family_id", ""),
            "gene_a": row.get("gene_a", ""),
            "gene_b": row.get("gene_b", ""),
            "Ks": "{:.8f}".format(value),
            "method": "CodeML_F3X4",
            "pair_class": "representative_A-B_collinear_homoeologue",
        })
    missing = set(AB_LABELS) - {str(row["comparison"]) for row in rows}
    if missing:
        raise ValueError("Panel a missing comparisons: {}".format(", ".join(sorted(missing))))
    return rows


def panel_b_rows(path: Path) -> list[dict[str, object]]:
    expected = {
        "{}_{}_{}".format(species, subgenome, reference)
        for species in AB_LABELS
        for subgenome in ("A", "B")
        for reference in ("PCAM", "PPE", "GTR")
    }
    rows: list[dict[str, object]] = []
    for row in read_tsv(path):
        species = row.get("species", "")
        subgenome = row.get("subgenome", "")
        reference = row.get("reference", "")
        label = "{}_{}_{}".format(species, subgenome, reference)
        if label not in expected:
            continue
        value = parse_ks(row.get("Ks", ""), path)
        rows.append({
            "comparison": label,
            "species": species,
            "subgenome": subgenome,
            "reference": reference,
            "pair_id": row.get("pair_id", ""),
            "query_gene": row.get("query_gene", ""),
            "reference_gene": row.get("ref_gene", ""),
            "block_id": row.get("block_id", ""),
            "Ks": "{:.8f}".format(value),
            "method": row.get("method", ""),
            "pair_class": "cross-species_collinear_orthologue",
        })
    observed = {str(row["comparison"]) for row in rows}
    if observed != expected:
        raise ValueError("Panel b missing comparisons: {}".format(", ".join(sorted(expected - observed))))
    methods = {str(row["method"]) for row in rows}
    if methods != {"YN"}:
        raise ValueError("Panel b should contain only YN Ks rows; found {}".format(sorted(methods)))
    return rows


def panel_c_rows(self_root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for label in SUBGENOME_SELF + SINGLE_SELF:
        source = self_root / label / "audit" / "{}_SELF.all_pair_ks_audit.tsv".format(label)
        if not source.is_file() or source.stat().st_size == 0:
            raise ValueError("Missing self Ks audit: {}".format(source))
        if label in SUBGENOME_SELF:
            parameter_path = self_root / label / "run_parameters.tsv"
            if not parameter_path.is_file() or "analysis_scope\tindependent_subgenome_self_synteny" not in parameter_path.read_text(encoding="utf-8", errors="replace"):
                raise ValueError("Panel c requires an independently run self-synteny analysis for {}".format(label))
        statuses: dict[str, int] = defaultdict(int)
        for row in read_tsv(source):
            status = row.get("ks_status", "unknown")
            statuses[status] += 1
            if status != "retained_0_lt_Ks_lt_5.0":
                continue
            value = parse_ks(row.get("dS", ""), source)
            rows.append({
                "comparison": label,
                "family_id": row.get("family_id", ""),
                "gene_1": row.get("copy_1", row.get("gene1", "")),
                "gene_2": row.get("copy_2", row.get("gene2", "")),
                "Ks": "{:.8f}".format(value),
                "method": "CodeML_F3X4",
                "pair_class": "independently_run_subgenome_cross-chromosome_self-syntenic_paralogue"
                if label in SUBGENOME_SELF else "genome-native_cross-chromosome_self-syntenic_paralogue",
            })
        for status, count in sorted(statuses.items()):
            audit_rows.append({
                "comparison": label,
                "source_dataset": label + "_SELF",
                "ks_status": status,
                "pairs": count,
            })

    expected = set(SUBGENOME_SELF) | set(SINGLE_SELF)
    observed = {str(row["comparison"]) for row in rows}
    if observed != expected:
        raise ValueError("Panel c missing retained comparisons: {}".format(", ".join(sorted(expected - observed))))
    return rows, audit_rows


def summaries(panel: str, rows: list[dict[str, object]]) -> list[dict[str, object]]:
    values: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        values[str(row["comparison"])].append(float(row["Ks"]))
    output: list[dict[str, object]] = []
    for comparison in sorted(values):
        subset = values[comparison]
        output.append({
            "panel": panel,
            "comparison": comparison,
            "n": len(subset),
            "Ks_median": "{:.8f}".format(statistics.median(subset)),
            "Ks_mean": "{:.8f}".format(statistics.mean(subset)),
            "filter": "0<Ks<5",
        })
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ab", type=Path, required=True, help="Formal CodeML A-B homoeologue Ks table")
    parser.add_argument("--reference", type=Path, required=True, help="Completed 18-comparison YN table")
    parser.add_argument("--self-root", type=Path, required=True, help="Completed CodeML self-synteny root")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    a_rows = panel_a_rows(args.ab)
    b_rows = panel_b_rows(args.reference)
    c_rows, c_audit = panel_c_rows(args.self_root)

    write_tsv(args.outdir / "panel_a_ab_homoeologue_ks.tsv",
              ["comparison", "family_id", "gene_a", "gene_b", "Ks", "method", "pair_class"], a_rows)
    write_tsv(args.outdir / "panel_b_subgenome_reference_ks.tsv",
              ["comparison", "species", "subgenome", "reference", "pair_id", "query_gene", "reference_gene", "block_id", "Ks", "method", "pair_class"], b_rows)
    write_tsv(args.outdir / "panel_c_subgenome_self_ks.tsv",
              ["comparison", "family_id", "gene_1", "gene_2", "Ks", "method", "pair_class"], c_rows)
    write_tsv(args.outdir / "panel_c_subgenome_self_filter_audit.tsv",
              ["comparison", "source_dataset", "ks_status", "pairs"], c_audit)
    summary_rows = summaries("a", a_rows) + summaries("b", b_rows) + summaries("c", c_rows)
    write_tsv(args.outdir / "input_summary.tsv",
              ["panel", "comparison", "n", "Ks_median", "Ks_mean", "filter"], summary_rows)
    (args.outdir / "method_contract.txt").write_text(
        "Panel a: existing representative PZI A1-B1, PSP A1-B1 and PPA A-B i-ADHoRe collinear homoeologues; CodeML F3X4 Ks.\n"
        "Panel b: completed independent A/B versus PCAM, PPE and GTR JCVI-collinear orthologues; KaKs_Calculator YN Ks.\n"
        "Panel c: independently run PZI_A, PZI_B, PSP_A, PSP_B, PPA_A and PPA_B CodeML self-synteny analyses; PCAM, PPE and GTR remain native self-synteny sets.\n"
        "All input rows passed the originating workflows' 0<Ks<5 filter. No Ks estimate or pair selection is recalculated during this preparation step.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Three-panel piecewise Ks inputs prepared.\n", encoding="ascii")
    print("[OK] Prepared panel a={}, panel b={}, panel c={} Ks rows in {}".format(
        len(a_rows), len(b_rows), len(c_rows), args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
