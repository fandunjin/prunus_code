#!/usr/bin/env python3
"""Audit the PPE self-synteny input and its low-Ks CodeML signal.

This is a non-plotting audit. It separates model redundancy, exact peptide
duplicates, cross-chromosome anchor structure, and the retained Ks tail so
that a low-Ks mode is not interpreted as a WGD without checking the inputs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Missing header: {path}")
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_bed(path: Path) -> tuple[dict[str, dict[str, str]], dict[tuple[str, str, str, str], list[str]]]:
    by_id: dict[str, dict[str, str]] = {}
    by_coord: defaultdict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise RuntimeError(f"Malformed BED line {line_number}: {path}")
            record = {
                "seqid": fields[0],
                "start": fields[1],
                "end": fields[2],
                "id": fields[3],
                "strand": fields[5],
            }
            if record["id"] in by_id:
                raise RuntimeError(f"Duplicate BED identifier {record['id']}: {path}")
            by_id[record["id"]] = record
            by_coord[(fields[0], fields[1], fields[2], fields[5])].append(record["id"])
    return by_id, dict(by_coord)


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in sequences:
                    raise RuntimeError(f"Duplicate FASTA identifier {current}: {path}")
                sequences[current] = []
            elif current is not None:
                sequences[current].append(line)
    return {key: "".join(value) for key, value in sequences.items()}


def finite_ks(row: dict[str, str]) -> float | None:
    try:
        value = float(row.get("dS", ""))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def median(values: list[float]) -> str:
    if not values:
        return "NA"
    values = sorted(values)
    middle = len(values) // 2
    if len(values) % 2:
        return f"{values[middle]:.8g}"
    return f"{((values[middle - 1] + values[middle]) / 2):.8g}"


def kde_mode(values: list[float], lower: float = 0.02, upper: float = 3.0,
             bandwidth: float = 0.035) -> tuple[str, str]:
    """Return the fixed-bandwidth self-panel KDE mode on the original Ks scale."""
    values = [value for value in values if lower <= value <= upper]
    if not values:
        return "NA", "NA"
    grid_count = int(round((upper - lower) / 0.002)) + 1
    density_factor = 1.0 / (math.sqrt(2.0 * math.pi) * bandwidth * len(values))
    best_x = lower
    best_density = -1.0
    for index in range(grid_count):
        x_value = lower + index * (upper - lower) / (grid_count - 1)
        density = density_factor * sum(
            math.exp(-0.5 * ((x_value - value) / bandwidth) ** 2) for value in values
        )
        if density > best_density:
            best_x, best_density = x_value, density
    return f"{best_x:.8g}", f"{best_density:.8g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--pep", type=Path, required=True)
    parser.add_argument("--ks", type=Path, required=True,
                        help="Retained 0<Ks<5 CodeML audit table")
    parser.add_argument("--all-pair-ks", type=Path,
                        help="Optional all-pair CodeML audit table including exclusions")
    parser.add_argument("--pairs", type=Path, required=True,
                        help="i-ADHoRe self_collinear_pairs.tsv")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.bed, args.pep, args.ks, args.pairs):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing input: {path}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    by_id, by_coord = read_bed(args.bed)
    peptides = read_fasta(args.pep)
    ks_rows = read_tsv(args.ks)
    pair_rows = read_tsv(args.pairs)

    duplicate_ids = {gene_id for ids in by_coord.values() if len(ids) > 1 for gene_id in ids}
    duplicate_coord_groups = [ids for ids in by_coord.values() if len(ids) > 1]
    sequence_groups: defaultdict[str, list[str]] = defaultdict(list)
    for gene_id, sequence in peptides.items():
        sequence_groups[sequence].append(gene_id)
    duplicate_sequence_groups = [ids for ids in sequence_groups.values() if len(ids) > 1]
    duplicate_sequence_ids = {gene_id for ids in duplicate_sequence_groups for gene_id in ids}

    model_rows = [{
        "metric": "bed_records", "value": len(by_id), "detail": "all BED models",
    }, {
        "metric": "unique_coordinate_groups", "value": len(by_coord),
        "detail": "seqid/start/end/strand groups",
    }, {
        "metric": "duplicate_coordinate_groups", "value": len(duplicate_coord_groups),
        "detail": "groups containing >1 model ID",
    }, {
        "metric": "models_in_duplicate_coordinate_groups", "value": len(duplicate_ids),
        "detail": "model IDs in duplicate coordinate groups",
    }, {
        "metric": "peptide_records", "value": len(peptides), "detail": "FASTA records",
    }, {
        "metric": "exact_duplicate_peptide_groups", "value": len(duplicate_sequence_groups),
        "detail": "groups with identical peptide sequence",
    }, {
        "metric": "records_in_exact_duplicate_peptide_groups", "value": len(duplicate_sequence_ids),
        "detail": "FASTA records in exact duplicate groups",
    }]
    write_tsv(args.outdir / "ppe_model_sequence_audit.tsv", model_rows, ["metric", "value", "detail"])

    thresholds = (0.02, 0.05, 0.10, 0.20, 0.50, 1.00)
    threshold_rows: list[dict[str, object]] = []
    low_examples: list[dict[str, object]] = []
    retained_records: list[tuple[dict[str, str], float, dict[str, str] | None,
                                 dict[str, str] | None, bool, bool]] = []
    for row in ks_rows:
        value = finite_ks(row)
        if value is None or not (0.0 < value < 5.0):
            continue
        first = by_id.get(row.get("gene_1", ""))
        second = by_id.get(row.get("gene_2", ""))
        duplicate_model = bool(
            row.get("gene_1", "") in duplicate_ids
            or row.get("gene_2", "") in duplicate_ids
        )
        exact_peptide = bool(
            row.get("gene_1", "") in duplicate_sequence_ids
            or row.get("gene_2", "") in duplicate_sequence_ids
        )
        retained_records.append((row, value, first, second, duplicate_model, exact_peptide))
    for threshold in thresholds:
        subset = []
        for row, value, first, second, duplicate_model, exact_peptide in retained_records:
            if not (value < threshold):
                continue
            subset.append((row, value, first, second, duplicate_model, exact_peptide))
        duplicate_count = sum(item[4] for item in subset)
        exact_count = sum(item[5] for item in subset)
        identities = [float(item[0]["alignmentidentity"]) for item in subset
                      if item[0].get("alignmentidentity", "")]
        coverages = [float(item[0]["alignmentcoverage"]) for item in subset
                     if item[0].get("alignmentcoverage", "")]
        threshold_rows.append({
            "threshold_Ks": f"{threshold:g}",
            "pairs": len(subset),
            "duplicate_coordinate_model_pairs": duplicate_count,
            "duplicate_coordinate_fraction": f"{(duplicate_count / len(subset)) if subset else 0:.8g}",
            "exact_duplicate_peptide_pairs": exact_count,
            "exact_duplicate_peptide_fraction": f"{(exact_count / len(subset)) if subset else 0:.8g}",
            "median_alignment_identity": median(identities),
            "median_alignment_coverage": median(coverages),
        })
        if threshold == 0.10:
            for row, value, first, second, duplicate_model, exact_peptide in sorted(subset, key=lambda item: item[1])[:50]:
                low_examples.append({
                    "Ks": f"{value:.8g}", "family_id": row.get("family_id", ""),
                    "gene_1": row.get("gene_1", ""), "gene_2": row.get("gene_2", ""),
                    "seqid_1": first["seqid"] if first else "NA",
                    "seqid_2": second["seqid"] if second else "NA",
                    "alignmentidentity": row.get("alignmentidentity", ""),
                    "alignmentcoverage": row.get("alignmentcoverage", ""),
                    "duplicate_coordinate_model": str(duplicate_model),
                    "exact_duplicate_peptide": str(exact_peptide),
                })
    write_tsv(args.outdir / "ppe_ks_threshold_audit.tsv", threshold_rows, list(threshold_rows[0]))
    write_tsv(args.outdir / "ppe_low_ks_examples.tsv", low_examples, list(low_examples[0]))

    sensitivity_specs = (
        ("original", lambda item: True),
        ("drop_duplicate_coordinate_models", lambda item: not item[4]),
        ("drop_exact_duplicate_peptide_models", lambda item: not item[5]),
        ("drop_either_model_risk", lambda item: not item[4] and not item[5]),
    )
    sensitivity_rows: list[dict[str, object]] = []
    for label, keep in sensitivity_specs:
        values = [value for item in retained_records if keep(item)
                  for value in [item[1]] if 0.02 <= value <= 3.0]
        mode_ks, mode_density = kde_mode(values)
        sensitivity_rows.append({
            "filter": label,
            "pairs_0.02_to_3": len(values),
            "pairs_Ks_lt_0.1": sum(value < 0.1 for value in values),
            "mode_Ks": mode_ks,
            "mode_density": mode_density,
        })
    write_tsv(args.outdir / "ppe_filter_sensitivity.tsv", sensitivity_rows, list(sensitivity_rows[0]))

    pair_audit = [{
        "metric": "iADHoRe_pair_records", "value": len(pair_rows),
        "detail": "self_collinear_pairs.tsv records",
    }, {
        "metric": "unique_pair_ids", "value": len({row.get("gene_1", "") + "__" + row.get("gene_2", "") for row in pair_rows}),
        "detail": "deduplicated gene-pair identifiers",
    }, {
        "metric": "cross_chromosome_records", "value": sum(row.get("seqid_1") != row.get("seqid_2") for row in pair_rows),
        "detail": "records retained by the self-synteny extractor",
    }, {
        "metric": "distinct_multiplicons", "value": len({row.get("multiplicon", "") for row in pair_rows}),
        "detail": "multiplicon identifiers",
    }]
    write_tsv(args.outdir / "ppe_pair_audit.tsv", pair_audit, ["metric", "value", "detail"])

    retained_values = [value for row in ks_rows if (value := finite_ks(row)) is not None and 0.0 < value < 5.0]
    all_pair_summary = ""
    if args.all_pair_ks and args.all_pair_ks.is_file():
        all_rows = read_tsv(args.all_pair_ks)
        saturated = sum((value := finite_ks(row)) is not None and value >= 5.0 for row in all_rows)
        nonpositive = sum((value := finite_ks(row)) is not None and value <= 0.0 for row in all_rows)
        all_pair_summary = f"\nAll-pair audit rows={len(all_rows)}; Ks>=5={saturated}; non-positive finite Ks={nonpositive}."
    report = (
        "PPE self-synteny Ks audit\n"
        "========================\n"
        f"BED models={len(by_id)}; unique coordinate groups={len(by_coord)}; "
        f"duplicate coordinate groups={len(duplicate_coord_groups)}; models in duplicate groups={len(duplicate_ids)}.\n"
        f"Peptide records={len(peptides)}; exact duplicate peptide groups={len(duplicate_sequence_groups)}; "
        f"records in exact duplicate groups={len(duplicate_sequence_ids)}.\n"
        f"Retained CodeML rows={len(ks_rows)}; 0<Ks<5 values={len(retained_values)}; "
        f"median Ks={median(retained_values)}.\n"
        f"Ks<0.1 is summarized in ppe_ks_threshold_audit.tsv; low-Ks examples are limited to the first 50 sorted pairs."
        f"{all_pair_summary}\n"
        "Interpretation: duplicate models and exact duplicate peptides are input risks, not evidence of WGD. "
        "A low-Ks mode should only be considered after this audit and an independent de-redundancy sensitivity run.\n"
    )
    (args.outdir / "ppe_audit_report.txt").write_text(report, encoding="ascii")
    print(f"[OK] PPE self-Ks audit written: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
