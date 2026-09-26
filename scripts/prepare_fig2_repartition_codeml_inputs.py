#!/usr/bin/env python3
"""Prepare one two-gene CodeML family set per corrected Figure 2 comparison.

The script intentionally keeps every phased comparison separate.  PZI and PSP
haplotypes are pooled only after pair-level Ks estimation, never by combining
haplotypes before CodeML.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


GROUPS = {
    "ab": "ab_pair_plan.tsv",
    "external": "external_pair_plan.tsv",
    "self": "self_pair_plan.tsv",
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


def read_fasta(path):
    records = {}
    record_id = None
    sequence = []

    def commit():
        if record_id is None:
            return
        if record_id in records:
            raise RuntimeError("Duplicate FASTA identifier in {}: {}".format(path, record_id))
        value = "".join(sequence).replace(" ", "").upper()
        if not value:
            raise RuntimeError("Empty FASTA sequence in {}: {}".format(path, record_id))
        records[record_id] = value

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                record_id = line[1:].split(None, 1)[0]
                sequence = []
            else:
                if record_id is None:
                    raise RuntimeError("FASTA sequence precedes a header in {}".format(path))
                sequence.append(line)
    commit()
    if not records:
        raise RuntimeError("No FASTA records in {}".format(path))
    return records


def write_fasta(path, records):
    with path.open("w", encoding="ascii", newline="") as handle:
        for gene in sorted(records):
            handle.write(">{}\n".format(gene))
            sequence = records[gene]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError("Missing input: {}".format(path))


def jcvi_pair_records(base, plan):
    pair_table = base / "01_jcvi_pairs" / "jcvi_syntenic_pairs.tsv"
    require_file(pair_table)
    by_comparison = defaultdict(list)
    for row in read_tsv(pair_table):
        by_comparison[row.get("comparison", "")].append(row)

    result = []
    for item in plan:
        comparison = item.get("comparison", "")
        query = item.get("query_label", "")
        reference = item.get("ref_label", "")
        if not comparison or not query or not reference:
            raise RuntimeError("Invalid pair plan row: {}".format(item))
        source_rows = by_comparison.get(comparison, [])
        if not source_rows:
            raise RuntimeError("No corrected JCVI anchors for {}".format(comparison))
        records = []
        seen_pairs = set()
        for index, row in enumerate(source_rows, start=1):
            if row.get("query_label") != query or row.get("ref_label") != reference:
                raise RuntimeError("JCVI label mismatch in {}".format(comparison))
            gene_1, gene_2 = row.get("query_gene", ""), row.get("ref_gene", "")
            if not gene_1 or not gene_2 or gene_1 == gene_2:
                raise RuntimeError("Invalid JCVI anchor in {}".format(comparison))
            pair_key = (gene_1, gene_2)
            if pair_key in seen_pairs:
                raise RuntimeError("Duplicate directed JCVI anchor in {}: {}".format(comparison, pair_key))
            seen_pairs.add(pair_key)
            records.append({
                "family_id": "{}_{:07d}".format(comparison, index),
                "gene_1": gene_1,
                "gene_2": gene_2,
                "source_pair_id": "{}:{}".format(comparison, index),
                "source_block_id": row.get("block_id", ""),
                "query_label": query,
                "ref_label": reference,
                "query_seqid": row.get("query_seqid", ""),
                "ref_seqid": row.get("ref_seqid", ""),
            })
        result.append((item, records, pair_table))
    return result


def self_pair_records(base, plan):
    result = []
    for item in plan:
        label = item.get("label", "")
        if not label:
            raise RuntimeError("Self plan contains an empty label")
        source = base / "02_iadhore_self" / label / "03_pairs" / "{}.families.tsv".format(label + "_SELF")
        require_file(source)
        source_rows = read_tsv(source)
        records = []
        seen_pairs = set()
        for index, row in enumerate(source_rows, start=1):
            source_id = row.get("family_id", "")
            gene_1, gene_2 = row.get("copy_1", ""), row.get("copy_2", "")
            if not source_id or not gene_1 or not gene_2 or gene_1 == gene_2:
                raise RuntimeError("Invalid self-synteny pair in {}".format(source))
            pair_key = tuple(sorted((gene_1, gene_2)))
            if pair_key in seen_pairs:
                raise RuntimeError("Duplicate self-synteny pair in {}: {}".format(label, pair_key))
            seen_pairs.add(pair_key)
            records.append({
                "family_id": "{}_self_{:07d}".format(label, index),
                "gene_1": gene_1,
                "gene_2": gene_2,
                "source_pair_id": source_id,
                "source_block_id": "",
                "query_label": label,
                "ref_label": label,
                "query_seqid": "",
                "ref_seqid": "",
            })
        if not records:
            raise RuntimeError("No self-synteny families in {}".format(source))
        comparison = "{}_self".format(label)
        planned = {
            "comparison": comparison,
            "query_label": label,
            "ref_label": label,
            "comparison_type": "self_synteny",
            "panel": item.get("panel", ""),
            "curve": item.get("curve", ""),
        }
        result.append((planned, records, source))
    return result


def write_comparison(outdir, planned, records, source, inputs_root, group):
    comparison = planned["comparison"]
    comparison_out = outdir / comparison
    comparison_out.mkdir(parents=True, exist_ok=False)
    labels = [planned["query_label"]]
    if planned["ref_label"] != planned["query_label"]:
        labels.append(planned["ref_label"])

    needed = defaultdict(set)
    for row in records:
        needed[row["query_label"]].add(row["gene_1"])
        needed[row["ref_label"]].add(row["gene_2"])

    selected = {}
    provenance = {}
    for label in labels:
        cds = inputs_root / label / (label + ".cds")
        require_file(cds)
        all_records = read_fasta(cds)
        missing = sorted(needed[label] - set(all_records))
        if missing:
            raise RuntimeError("{} CDS IDs absent from {}: {}".format(
                len(missing), cds, ",".join(missing[:8])))
        for gene in needed[label]:
            if gene in selected:
                raise RuntimeError("CDS identifier occurs in multiple inputs: {}".format(gene))
            selected[gene] = all_records[gene]
            provenance[gene] = label

    # WGD receives a strict three-column family table. The full per-pair
    # provenance remains in pair_manifest.tsv below, rather than leaking
    # additional fields into the WGD input.
    write_tsv(
        comparison_out / "families.tsv",
        ["family_id", "gene_1", "gene_2"],
        [{"family_id": row["family_id"], "gene_1": row["gene_1"], "gene_2": row["gene_2"]}
         for row in records],
    )
    write_tsv(
        comparison_out / "pair_manifest.tsv",
        [
            "family_id", "gene_1", "gene_2", "source_pair_id", "source_block_id",
            "query_label", "ref_label", "query_seqid", "ref_seqid",
        ],
        records,
    )
    write_fasta(comparison_out / "combined.cds.fa", selected)
    write_tsv(
        comparison_out / "cds_manifest.tsv",
        ["gene", "source_label", "cds_length"],
        [{"gene": gene, "source_label": provenance[gene], "cds_length": len(selected[gene])} for gene in sorted(selected)],
    )
    write_tsv(
        comparison_out / "preparation_summary.tsv",
        ["comparison", "group", "source", "families", "unique_cds", "query_label", "ref_label"],
        [{
            "comparison": comparison,
            "group": group,
            "source": str(source),
            "families": len(records),
            "unique_cds": len(selected),
            "query_label": planned["query_label"],
            "ref_label": planned["ref_label"],
        }],
    )
    return {
        "comparison": comparison,
        "group": group,
        "query_label": planned["query_label"],
        "ref_label": planned["ref_label"],
        "comparison_type": planned.get("comparison_type", ""),
        "panel": planned.get("panel", ""),
        "curve": planned.get("curve", ""),
        "families": len(records),
        "unique_cds": len(selected),
        "source": str(source),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--group", choices=sorted(GROUPS), required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    base = args.base.resolve()
    inputs = base / "inputs"
    plan_path = base / "plans" / GROUPS[args.group]
    require_file(plan_path)
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    plan = read_tsv(plan_path)
    if args.group == "self":
        prepared = self_pair_records(base, plan)
    else:
        prepared = jcvi_pair_records(base, plan)

    manifest = []
    for planned, records, source in prepared:
        manifest.append(write_comparison(args.outdir, planned, records, source, inputs, args.group))
    if not manifest:
        raise SystemExit("[ERROR] No comparisons were prepared")
    write_tsv(
        args.outdir / "comparisons.tsv",
        ["comparison", "group", "query_label", "ref_label", "comparison_type", "panel", "curve", "families", "unique_cds", "source"],
        manifest,
    )
    (args.outdir / "method_contract.txt").write_text(
        "Input pairs are phase-specific corrected JCVI anchors (A-B/external) or i-ADHoRe cross-chromosome anchorpoints (self).\n"
        "Each output family has exactly two genes and is estimated independently by CodeML.\n"
        "PZI and PSP haplotypes remain separate until the plotting aggregation stage.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("CodeML family inputs prepared.\n", encoding="ascii")
    print("[OK] Prepared {} {} comparisons: {}".format(len(manifest), args.group, args.outdir))


if __name__ == "__main__":
    main()
