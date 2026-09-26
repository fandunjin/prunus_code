#!/usr/bin/env python3
"""Extract audit-ready structural-junction candidates from Fig. 4 WGDI maps.

This script labels only neutral evidence classes: a source-unit transition,
an orientation transition, or both. Mechanistic labels (RCT, EEJ, inversion)
are intentionally left blank until a candidate is corroborated by the
cross-species and outgroup comparisons specified in the curation template.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Gene:
    chromosome: str
    order: int
    start: int
    end: int
    ident: str

    @property
    def center(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass(frozen=True)
class Block:
    identifier: str
    source_unit: str
    target_chromosome: str
    source_start: int
    source_end: int
    target_start: int
    target_end: int
    length: int
    density1: float
    density2: float
    pvalue: float

    @property
    def target_low(self) -> int:
        return min(self.target_start, self.target_end)

    @property
    def target_high(self) -> int:
        return max(self.target_start, self.target_end)

    @property
    def source_low(self) -> int:
        return min(self.source_start, self.source_end)

    @property
    def source_high(self) -> int:
        return max(self.source_start, self.source_end)

    @property
    def orientation(self) -> str:
        return "+" if (self.source_end - self.source_start) * (self.target_end - self.target_start) >= 0 else "-"

    @property
    def rank(self) -> tuple[float, ...]:
        return (self.length, self.density1 + self.density2, -self.pvalue)


def integer(value: str) -> int:
    return int(float(value))


def number(value: str, default: float = 0.0) -> float:
    try:
        parsed = float(value)
        return parsed if math.isfinite(parsed) else default
    except (TypeError, ValueError):
        return default


def read_gff(path: Path) -> tuple[dict[str, dict[int, Gene]], dict[str, int]]:
    by_chromosome: dict[str, dict[int, Gene]] = defaultdict(dict)
    lengths: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                start, end, order = integer(fields[2]), integer(fields[3]), integer(fields[5])
            except ValueError:
                continue
            start, end = sorted((start, end))
            gene = Gene(fields[0], order, start, end, fields[1])
            by_chromosome[fields[0]][order] = gene
            lengths[fields[0]] = max(lengths[fields[0]], end)
    if not by_chromosome:
        raise ValueError(f"No usable GFF rows: {path}")
    return by_chromosome, lengths


def read_blocks(path: Path, allowed_units: set[str], min_length: int) -> list[Block]:
    blocks: list[Block] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for index, row in enumerate(csv.DictReader(handle), 1):
            if row.get("chr1", "") not in allowed_units:
                continue
            try:
                length = integer(row.get("length", "0"))
                start1, end1 = integer(row["start1"]), integer(row["end1"])
                start2, end2 = integer(row["start2"]), integer(row["end2"])
            except (KeyError, ValueError):
                continue
            if length < min_length:
                continue
            blocks.append(Block(
                identifier=row.get("id", str(index)), source_unit=row["chr1"], target_chromosome=row.get("chr2", ""),
                source_start=start1, source_end=end1, target_start=start2, target_end=end2, length=length,
                density1=number(row.get("density1", "")), density2=number(row.get("density2", "")),
                pvalue=number(row.get("pvalue", ""), 1.0),
            ))
    return blocks


def overlaps(a_low: int, a_high: int, b_low: int, b_high: int) -> bool:
    return not (a_high < b_low or b_high < a_low)


def choose_primary(blocks: list[Block], max_target_overlap: float) -> list[Block]:
    """Greedily retain the best target-non-overlapping WGDI blocks.

    This removes nested alternative chains but does not discard non-overlapping
    uses of the same source unit, which may be biologically informative.
    """
    selected: list[Block] = []
    by_chromosome: dict[str, list[Block]] = defaultdict(list)
    for block in blocks:
        by_chromosome[block.target_chromosome].append(block)
    for chromosome, group in by_chromosome.items():
        accepted: list[Block] = []
        for block in sorted(group, key=lambda item: item.rank, reverse=True):
            reject = False
            for existing in accepted:
                overlap = max(0, min(block.target_high, existing.target_high) - max(block.target_low, existing.target_low) + 1)
                smaller = min(block.target_high - block.target_low + 1, existing.target_high - existing.target_low + 1)
                if smaller and overlap / smaller > max_target_overlap:
                    reject = True
                    break
            if not reject:
                accepted.append(block)
        selected.extend(sorted(accepted, key=lambda item: (item.target_low, item.target_high, item.source_unit)))
    return selected


def nearest_gene(order_map: dict[int, Gene], order: int, side: str) -> Gene | None:
    values = sorted(order_map)
    if side == "left":
        candidates = [value for value in values if value <= order]
        return order_map[candidates[-1]] if candidates else None
    candidates = [value for value in values if value >= order]
    return order_map[candidates[0]] if candidates else None


def target_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="0903 directory")
    parser.add_argument("--min-block-genes", type=int, default=10)
    parser.add_argument("--max-target-overlap", type=float, default=0.20)
    parser.add_argument("--max-gap-genes", type=int, default=100)
    args = parser.parse_args()
    if args.min_block_genes < 2 or not 0 <= args.max_target_overlap < 1 or args.max_gap_genes < 0:
        raise ValueError("Invalid block-selection parameters")

    root = args.root.resolve()
    work = root / "karyotype" / "fig4_breakpoint_candidates"
    reference_gff = work / "reference" / "prunus9_reference.gff"
    target_rows = target_metadata(work / "target_manifest.tsv")
    source_gff, _source_lengths = read_gff(reference_gff)
    allowed_units = set(source_gff)
    candidate_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []

    for metadata in target_rows:
        label = metadata["target_label"]
        blockinfo = work / "mapping" / label / f"prunus9_vs_{label}.blockinfo.csv"
        if not blockinfo.is_file() or blockinfo.stat().st_size == 0:
            continue
        target_gff, target_lengths = read_gff(root / metadata["gff"])
        all_blocks = read_blocks(blockinfo, allowed_units, args.min_block_genes)
        primary = choose_primary(all_blocks, args.max_target_overlap)
        for block in primary:
            selected_rows.append({
                "target_label": label, "species": metadata["species"], "subgenome": metadata["subgenome"],
                "target_chromosome": block.target_chromosome, "raw_block_id": block.identifier,
                "source_unit": block.source_unit, "source_start_order": block.source_start, "source_end_order": block.source_end,
                "target_start_order": block.target_start, "target_end_order": block.target_end, "orientation": block.orientation,
                "length_genes": block.length, "density_source": f"{block.density1:.6g}", "density_target": f"{block.density2:.6g}",
                "pvalue": f"{block.pvalue:.6g}",
            })
        by_target: dict[str, list[Block]] = defaultdict(list)
        for block in primary:
            by_target[block.target_chromosome].append(block)
        for chromosome, chain in by_target.items():
            if chromosome not in target_gff:
                continue
            chain.sort(key=lambda item: (item.target_low, item.target_high))
            order_map = target_gff[chromosome]
            chromosome_length = target_lengths.get(chromosome, 0)
            for left, right in zip(chain, chain[1:]):
                gap = right.target_low - left.target_high - 1
                if gap < 0 or gap > args.max_gap_genes:
                    continue
                source_change = left.source_unit != right.source_unit
                orientation_change = left.orientation != right.orientation
                if not source_change and not orientation_change:
                    continue
                left_gene = nearest_gene(order_map, left.target_high, "left")
                right_gene = nearest_gene(order_map, right.target_low, "right")
                if left_gene is None or right_gene is None:
                    continue
                low_bp, high_bp = sorted((int(left_gene.center), int(right_gene.center)))
                evidence_class = (
                    "source_and_orientation_transition" if source_change and orientation_change
                    else "source_transition" if source_change else "orientation_transition"
                )
                pair = "|".join(sorted((left.source_unit, right.source_unit)))
                candidate_id = f"{label}_{chromosome}_{left.identifier}_{right.identifier}"
                candidate_rows.append({
                    "candidate_id": candidate_id,
                    "evidence_class": evidence_class,
                    "target_label": label,
                    "species": metadata["species"],
                    "subgenome": metadata["subgenome"],
                    "target_chromosome": chromosome,
                    "target_chromosome_length_bp": chromosome_length,
                    "source_left_unit": left.source_unit,
                    "source_right_unit": right.source_unit,
                    "source_unit_pair": pair,
                    "left_orientation": left.orientation,
                    "right_orientation": right.orientation,
                    "left_block_id": left.identifier,
                    "right_block_id": right.identifier,
                    "left_block_gene_count": left.length,
                    "right_block_gene_count": right.length,
                    "interblock_gap_genes": gap,
                    "candidate_low_bp": low_bp,
                    "candidate_high_bp": high_bp,
                    "candidate_center_bp": round((low_bp + high_bp) / 2),
                    "candidate_bracket_width_bp": high_bp - low_bp,
                    "left_anchor_gene": left_gene.ident,
                    "left_anchor_start_bp": left_gene.start,
                    "left_anchor_end_bp": left_gene.end,
                    "right_anchor_gene": right_gene.ident,
                    "right_anchor_start_bp": right_gene.start,
                    "right_anchor_end_bp": right_gene.end,
                    "mechanism_status": "unassigned_pending_cross_species_and_outgroup_review",
                })

    selected_fields = [
        "target_label", "species", "subgenome", "target_chromosome", "raw_block_id", "source_unit",
        "source_start_order", "source_end_order", "target_start_order", "target_end_order", "orientation",
        "length_genes", "density_source", "density_target", "pvalue",
    ]
    candidate_fields = [
        "candidate_id", "evidence_class", "target_label", "species", "subgenome", "target_chromosome",
        "target_chromosome_length_bp", "source_left_unit", "source_right_unit", "source_unit_pair",
        "left_orientation", "right_orientation", "left_block_id", "right_block_id", "left_block_gene_count",
        "right_block_gene_count", "interblock_gap_genes", "candidate_low_bp", "candidate_high_bp",
        "candidate_center_bp", "candidate_bracket_width_bp", "left_anchor_gene", "left_anchor_start_bp",
        "left_anchor_end_bp", "right_anchor_gene", "right_anchor_start_bp", "right_anchor_end_bp", "mechanism_status",
    ]
    write_tsv(work / "selected_primary_blocks.tsv", selected_rows, selected_fields)
    write_tsv(work / "junction_candidates.tsv", candidate_rows, candidate_fields)

    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in candidate_rows:
        key = (str(row["evidence_class"]), str(row["source_unit_pair"]), str(row["left_orientation"]), str(row["right_orientation"]))
        grouped[key].append(row)
    summary_rows: list[dict[str, object]] = []
    for index, (key, rows) in enumerate(sorted(grouped.items()), 1):
        evidence_class, pair, left_orientation, right_orientation = key
        species = sorted({str(row["species"]) for row in rows})
        labels = sorted({str(row["target_label"]) for row in rows})
        subgenomes = sorted({str(row["subgenome"]) for row in rows})
        summary_rows.append({
            "cluster_id": f"J{index:04d}", "evidence_class": evidence_class, "source_unit_pair": pair,
            "left_orientation": left_orientation, "right_orientation": right_orientation,
            "candidate_count": len(rows), "species_count": len(species), "species": ",".join(species),
            "target_count": len(labels), "targets": ",".join(labels), "subgenomes": ",".join(subgenomes),
            "candidate_ids": ",".join(str(row["candidate_id"]) for row in rows),
            "curation_rule": "Assign RCT, EEJ or inversion only after checking reciprocal topology and all listed target/outgroup occurrences.",
        })
    summary_fields = [
        "cluster_id", "evidence_class", "source_unit_pair", "left_orientation", "right_orientation", "candidate_count",
        "species_count", "species", "target_count", "targets", "subgenomes", "candidate_ids", "curation_rule",
    ]
    write_tsv(work / "junction_candidate_clusters.tsv", summary_rows, summary_fields)

    template_rows = []
    for row in summary_rows:
        if int(row["species_count"]) < 2:
            continue
        template_rows.append({
            "event_id": "", "mechanism": "", "lineage_assignment": "", "cluster_id": row["cluster_id"],
            "candidate_ids": row["candidate_ids"], "expected_evidence": "reciprocal junctions for RCT; source-terminal joins for EEJ; paired orientation switches for inversion",
            "outgroup_check": "PCAM/PPE/GTR reviewed", "curation_status": "PENDING", "notes": "",
        })
    template_fields = ["event_id", "mechanism", "lineage_assignment", "cluster_id", "candidate_ids", "expected_evidence", "outgroup_check", "curation_status", "notes"]
    write_tsv(work / "breakpoint_event_curation_template.tsv", template_rows, template_fields)

    (work / "candidate_extraction_parameters.tsv").write_text(
        "parameter\tvalue\n"
        f"min_block_genes\t{args.min_block_genes}\n"
        f"max_target_overlap\t{args.max_target_overlap}\n"
        f"max_gap_genes\t{args.max_gap_genes}\n"
        "coordinate_definition\tCandidate interval is bracketed by the nearest flanking WGDI anchor-gene centres in the target GFF.\n"
        "mechanism_policy\tNeutral evidence classes only; no candidate is automatically labelled RCT, EEJ or inversion.\n",
        encoding="utf-8",
    )
    print(f"[OK] primary blocks: {len(selected_rows):,}")
    print(f"[OK] neutral junction candidates: {len(candidate_rows):,}")
    print(f"[OK] clusters: {len(summary_rows):,}")


if __name__ == "__main__":
    main()
