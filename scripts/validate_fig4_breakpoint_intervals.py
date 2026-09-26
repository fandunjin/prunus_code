#!/usr/bin/env python3
"""Validate AK1/AK3 breakpoint candidates from the 0903 WGDI evidence.

The WGDI extraction provides the two flanking anchor genes for each
source-unit/orientation transition.  This script converts those anchors into
an explicit inter-anchor interval on the extant chromosome, summarizes
reciprocal and cross-species support, and records the evidence boundary for
single-nucleotide validation.  It deliberately does not call a midpoint an
exact breakpoint: a nucleotide-level junction requires a FASTA sequence for
the source/reference chromosome and sequence alignment across the junction.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FOCAL_UNITS = {"C2B": "AK1", "GTR03": "AK3"}
UNIT_TO_AK = {
    "C2B": "AK1",
    "C5B": "AK2",
    "GTR03": "AK3",
    "GTR04": "AK4",
    "GTR05": "AK5",
    "GTR06": "AK6",
    "GTR07": "AK7",
    "GTR08": "AK8",
    "GTR09": "AK9",
}
RACEMOSE = {"PZI", "PSP", "PPA"}
EXTERNAL = {"PCAM", "PPE", "GTR", "RCH"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def as_int(row: dict[str, str], key: str) -> int:
    return int(float(row[key]))


def pair_key(row: dict[str, str]) -> tuple[str, str]:
    return (row["source_left_unit"], row["source_right_unit"])


def sorted_pair(row: dict[str, str]) -> str:
    return "|".join(sorted((row["source_left_unit"], row["source_right_unit"])))


def anchor_interval(row: dict[str, str]) -> tuple[int, int, int, str]:
    """Return the bases between the two anchor features.

    The anchors are first ordered by their physical centres.  The returned
    interval is the sequence after the upstream anchor and before the
    downstream anchor.  If the two features overlap, the interval is empty
    and the row is explicitly marked as overlapping rather than guessed.
    """
    left_start = as_int(row, "left_anchor_start_bp")
    left_end = as_int(row, "left_anchor_end_bp")
    right_start = as_int(row, "right_anchor_start_bp")
    right_end = as_int(row, "right_anchor_end_bp")
    left_center = (left_start + left_end) / 2
    right_center = (right_start + right_end) / 2
    if left_center <= right_center:
        upstream_end, downstream_start = left_end, right_start
    else:
        upstream_end, downstream_start = right_end, left_start
    low = min(upstream_end, downstream_start) + 1
    high = max(upstream_end, downstream_start) - 1
    if upstream_end < downstream_start:
        status = "gene_bounded_interval"
        gap = high - low + 1
    else:
        status = "overlapping_anchor_features"
        gap = 0
    return low, high, gap, status


def target_reciprocal(rows: list[dict[str, str]], row: dict[str, str]) -> bool:
    """Check whether both source directions occur on the same target."""
    target = (row["target_label"], row["target_chromosome"], sorted_pair(row))
    directions = {
        (r["source_left_unit"], r["source_right_unit"])
        for r in rows
        if (r["target_label"], r["target_chromosome"], sorted_pair(r)) == target
    }
    if row["source_left_unit"] == row["source_right_unit"]:
        return (row["left_orientation"], row["right_orientation"]) in directions and (
            (row["right_orientation"], row["left_orientation"]) in {
                (r["left_orientation"], r["right_orientation"])
                for r in rows
                if (r["target_label"], r["target_chromosome"], sorted_pair(r)) == target
            }
        )
    return (row["source_right_unit"], row["source_left_unit"]) in directions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903"))
    args = parser.parse_args()
    root = args.root.resolve()
    work = root / "karyotype" / "fig4_breakpoint_candidates"
    rows = read_tsv(work / "junction_candidates.tsv")
    candidate_to_cluster: dict[str, str] = {}
    clusters_path = work / "junction_candidate_clusters.tsv"
    if clusters_path.is_file():
        for cluster in read_tsv(clusters_path):
            for candidate_id in cluster.get("candidate_ids", "").split(","):
                if candidate_id:
                    candidate_to_cluster[candidate_id] = cluster["cluster_id"]
    focal = [
        row for row in rows
        if row["source_left_unit"] in FOCAL_UNITS or row["source_right_unit"] in FOCAL_UNITS
    ]
    if not focal:
        raise SystemExit("No AK1/AK3-associated candidates found in junction_candidates.tsv")

    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_pair[sorted_pair(row)].append(row)

    interval_rows: list[dict[str, object]] = []
    for row in focal:
        low, high, gap, interval_status = anchor_interval(row)
        pair = sorted_pair(row)
        species = sorted({r["species"] for r in by_pair[pair]})
        racemose = sorted(set(species) & RACEMOSE)
        external = sorted(set(species) & EXTERNAL)
        reciprocal = target_reciprocal(rows, row)
        if reciprocal and len(racemose) >= 2:
            validation = "cross_species_reciprocal_candidate"
        elif reciprocal:
            validation = "reciprocal_candidate"
        else:
            validation = "single_transition_candidate"
        if external:
            validation += ";outgroup_supported"
        left_ak = UNIT_TO_AK.get(row["source_left_unit"], row["source_left_unit"])
        right_ak = UNIT_TO_AK.get(row["source_right_unit"], row["source_right_unit"])
        canonical_ak_pair = "|".join(
            UNIT_TO_AK.get(unit, unit) for unit in pair.split("|")
        )
        pair_scope = (
            "same_focal_unit" if row["source_left_unit"] == row["source_right_unit"] and row["source_left_unit"] in FOCAL_UNITS
            else "focal_to_focal" if row["source_left_unit"] in FOCAL_UNITS and row["source_right_unit"] in FOCAL_UNITS
            else "focal_to_other"
        )
        interval_rows.append({
            "candidate_id": row["candidate_id"],
            "cluster_id": candidate_to_cluster.get(row["candidate_id"], "unclustered"),
            "candidate_class": row["evidence_class"],
            "source_unit_pair": pair,
            "source_left_unit": row["source_left_unit"],
            "source_right_unit": row["source_right_unit"],
            "source_left_AK": left_ak,
            "source_right_AK": right_ak,
            "source_AK_pair": f"{left_ak}|{right_ak}",
            "canonical_AK_pair": canonical_ak_pair,
            "pair_scope": pair_scope,
            "species": row["species"],
            "target_label": row["target_label"],
            "subgenome": row["subgenome"],
            "target_chromosome": row["target_chromosome"],
            "target_chromosome_length_bp": row["target_chromosome_length_bp"],
            "left_anchor_gene": row["left_anchor_gene"],
            "left_anchor_start_bp": row["left_anchor_start_bp"],
            "left_anchor_end_bp": row["left_anchor_end_bp"],
            "right_anchor_gene": row["right_anchor_gene"],
            "right_anchor_start_bp": row["right_anchor_start_bp"],
            "right_anchor_end_bp": row["right_anchor_end_bp"],
            "inter_anchor_low_bp": low,
            "inter_anchor_high_bp": high,
            "inter_anchor_gap_bp": gap,
            "coordinate_system": "1-based_inclusive",
            "wgdi_candidate_center_bp": row["candidate_center_bp"],
            "wgdi_candidate_width_bp": row["candidate_bracket_width_bp"],
            "left_orientation": row["left_orientation"],
            "right_orientation": row["right_orientation"],
            "reciprocal_target_support": "yes" if reciprocal else "no",
            "racemose_species_support": ",".join(racemose),
            "outgroup_support": ",".join(external),
            "validation_status": validation,
            "interval_status": interval_status,
            "exact_base_status": "not_validated_without_source_reference_FASTA",
            "exact_base_requirement": "sequence_alignment_across_target_and_source_junction",
        })

    interval_fields = list(interval_rows[0])
    write_tsv(work / "confirmed_breakpoint_intervals_AK1_AK3.tsv", interval_rows, interval_fields)
    reciprocal_rows = [
        row for row in interval_rows
        if row["pair_scope"] != "same_focal_unit"
        and row["reciprocal_target_support"] == "yes"
    ]
    write_tsv(
        work / "validated_reciprocal_breakpoint_intervals_AK1_AK3.tsv",
        reciprocal_rows,
        interval_fields,
    )

    # One compact summary per source-unit pair.  It is intended as the input
    # ledger for later gene/TE analyses, while the interval table preserves
    # every haplotype and anchor-gene observation.
    summary: list[dict[str, object]] = []
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in interval_rows:
        grouped[(str(row["source_unit_pair"]), str(row["target_label"]), str(row["subgenome"]), str(row["target_chromosome"]))].append(row)
    for key, group in sorted(grouped.items()):
        pair, target, subgenome, chrom = key
        widths = [int(r["inter_anchor_gap_bp"]) for r in group if str(r["interval_status"]) == "gene_bounded_interval"]
        summary.append({
            "source_unit_pair": pair,
            "source_AK_pair": "|".join(UNIT_TO_AK.get(x, x) for x in pair.split("|")),
            "target_label": target,
            "subgenome": subgenome,
            "target_chromosome": chrom,
            "candidate_count": len(group),
            "target_interval_min_bp": min(int(r["inter_anchor_low_bp"]) for r in group),
            "target_interval_max_bp": max(int(r["inter_anchor_high_bp"]) for r in group),
            "target_interval_span_bp": max(int(r["inter_anchor_high_bp"]) for r in group) - min(int(r["inter_anchor_low_bp"]) for r in group) + 1,
            "gene_bounded_gap_min_bp": min(widths) if widths else "NA",
            "gene_bounded_gap_max_bp": max(widths) if widths else "NA",
            "reciprocal_support": group[0]["reciprocal_target_support"],
            "racemose_species_support": group[0]["racemose_species_support"],
            "outgroup_support": group[0]["outgroup_support"],
            "status": "validated_gene_bounded_interval;single_base_pending_reference_FASTA",
        })
    summary_fields = list(summary[0]) if summary else []
    write_tsv(work / "breakpoint_validation_summary_AK1_AK3.tsv", summary, summary_fields)

    combo_rows: list[dict[str, object]] = []
    combo_groups: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in reciprocal_rows:
        combo_groups[(
            str(row["canonical_AK_pair"]),
            str(row["species"]),
            str(row["target_label"]),
            str(row["subgenome"]),
            str(row["target_chromosome"]),
        )].append(row)
    for key, group in sorted(combo_groups.items()):
        ak_pair, species, target_label, subgenome, chromosome = key
        combo_rows.append({
            "canonical_AK_pair": ak_pair,
            "species": species,
            "target_label": target_label,
            "subgenome": subgenome,
            "target_chromosome": chromosome,
            "candidate_count": len(group),
            "candidate_ids": ",".join(str(row["candidate_id"]) for row in group),
            "inter_anchor_low_bp": min(int(row["inter_anchor_low_bp"]) for row in group),
            "inter_anchor_high_bp": max(int(row["inter_anchor_high_bp"]) for row in group),
            "racemose_species_support": group[0]["racemose_species_support"],
            "outgroup_support": group[0]["outgroup_support"],
            "status": "reciprocal_gene_bounded_interval",
        })
    combo_fields = list(combo_rows[0]) if combo_rows else []
    write_tsv(
        work / "current_chromosome_combinations_AK1_AK3.tsv",
        combo_rows,
        combo_fields,
    )

    methods = (
        "AK1=C2B and AK3=GTR03 were restored from the 0903 reference-unit mapping.\n"
        "Coordinates in confirmed_breakpoint_intervals_AK1_AK3.tsv are the inter-anchor intervals between the two flanking WGDI anchor genes.\n"
        "wgdi_candidate_center_bp is a midpoint of the original anchor-gene centres and is not an exact nucleotide breakpoint.\n"
        "A single-base breakpoint requires a sequence-level alignment to the corresponding source/reference chromosome FASTA; 0903 contains GFF/peptide inputs for GTR and outgroups but no corresponding chromosome FASTA.\n"
        "Rows remain evidence-labelled; RCT, EEJ and inversion assignments require reciprocal topology and outgroup curation.\n"
    )
    (work / "breakpoint_validation_methods_AK1_AK3.txt").write_text(methods, encoding="utf-8")
    print(f"[OK] AK1/AK3 interval rows: {len(interval_rows):,}")
    print(f"[OK] AK1/AK3 reciprocal non-self rows: {len(reciprocal_rows):,}")
    print(f"[OK] AK1/AK3 reciprocal chromosome combinations: {len(combo_rows):,}")
    print(f"[OK] AK1/AK3 target summaries: {len(summary):,}")


if __name__ == "__main__":
    main()
