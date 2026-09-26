#!/usr/bin/env python3
"""Extract conservative cross-chromosome self-syntenic anchor pairs from i-ADHoRe."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("No header in {}".format(path))
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--iadhore-input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    input_dir = args.iadhore_input
    multi_path = input_dir / "iadhore_output" / "multiplicons.txt"
    anchor_path = input_dir / "iadhore_output" / "anchorpoints.txt"
    list_path = input_dir / "list_manifest.tsv"
    for path in (multi_path, anchor_path, list_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit("[ERROR] Missing i-ADHoRe output/input: {}".format(path))
    args.outdir.mkdir(parents=True, exist_ok=True)

    list_to_seqid = {row["list_name"]: row["seqid"] for row in read_tsv(list_path)}
    cross_multiplicons: dict[str, tuple[str, str]] = {}
    for row in read_tsv(multi_path):
        if row.get("genome_x") != args.label or row.get("genome_y") != args.label:
            continue
        left = row.get("list_x", "")
        right = row.get("list_y", "")
        if not left or not right or left == right:
            continue
        if left not in list_to_seqid or right not in list_to_seqid:
            raise SystemExit("[ERROR] i-ADHoRe multiplicon refers to unrecognised list")
        cross_multiplicons[row["id"]] = (left, right)
    if not cross_multiplicons:
        raise SystemExit("[ERROR] No cross-chromosome self-syntenic multiplicons for {}".format(args.label))

    deduplicated: dict[tuple[str, str], dict[str, object]] = {}
    real_anchorpoints = 0
    for row in read_tsv(anchor_path):
        multiplicon = row.get("multiplicon", "")
        if multiplicon not in cross_multiplicons:
            continue
        if str(row.get("is_real_anchorpoint", "")).lower() in {"0", "false", "no"}:
            continue
        gene_x, gene_y = row.get("gene_x", ""), row.get("gene_y", "")
        if not gene_x or not gene_y or gene_x == gene_y:
            continue
        real_anchorpoints += 1
        left_list, right_list = cross_multiplicons[multiplicon]
        if gene_x < gene_y:
            key = (gene_x, gene_y)
            seq_1, seq_2 = list_to_seqid[left_list], list_to_seqid[right_list]
        else:
            key = (gene_y, gene_x)
            seq_1, seq_2 = list_to_seqid[right_list], list_to_seqid[left_list]
        previous = deduplicated.get(key)
        if previous is None:
            deduplicated[key] = {
                "label": args.label, "gene_1": key[0], "gene_2": key[1],
                "seqid_1": seq_1, "seqid_2": seq_2, "multiplicon": multiplicon,
                "duplicate_anchor_records": 1,
            }
        else:
            previous["duplicate_anchor_records"] = int(previous["duplicate_anchor_records"]) + 1
    if not deduplicated:
        raise SystemExit("[ERROR] No real cross-chromosome anchor pairs retained for {}".format(args.label))

    pair_rows = []
    family_rows = []
    for index, row in enumerate(sorted(deduplicated.values(), key=lambda item: (str(item["gene_1"]), str(item["gene_2"]))), 1):
        family_id = "{}_SELF_{:06d}".format(args.label, index)
        row["family_id"] = family_id
        pair_rows.append(row)
        family_rows.append({"family_id": family_id, "copy_1": row["gene_1"], "copy_2": row["gene_2"]})
    write_tsv(args.outdir / "self_collinear_pairs.tsv", [
        "family_id", "label", "gene_1", "gene_2", "seqid_1", "seqid_2", "multiplicon", "duplicate_anchor_records"], pair_rows)
    write_tsv(args.outdir / "{}_SELF.families.tsv".format(args.label), ["family_id", "copy_1", "copy_2"], family_rows)
    summary = [{"label": args.label, "cross_chromosome_multiplicons": len(cross_multiplicons),
                "real_cross_chromosome_anchor_records": real_anchorpoints,
                "unique_self_collinear_pairs": len(pair_rows)}]
    write_tsv(args.outdir / "self_collinear_pair_summary.tsv", list(summary[0]), summary)
    (args.outdir / "method_contract.txt").write_text(
        "Only real i-ADHoRe anchorpoints in within-genome multiplicons spanning two different chromosomes are retained.\n"
        "Same-list anchors, self-hits and duplicate pair records are excluded or deduplicated before pairwise Ks calculation.\n", encoding="ascii")
    (args.outdir / "COMPLETE").write_text("Cross-chromosome self-syntenic pairs extracted.\n", encoding="ascii")
    print("[OK] {}: {} unique pairs from {} cross-chromosome multiplicons".format(
        args.label, len(pair_rows), len(cross_multiplicons)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
