#!/usr/bin/env python3
"""Summarize A/B syntenic anchor states against a low-copy reference.

This script intentionally distinguishes an absent JCVI anchor from a verified
gene loss. It creates a candidate table for later sequence-level validation;
it never labels an unanchored locus as a deletion.
"""

from __future__ import print_function

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


STATE_ORDER = [
    "both_anchored",
    "A_only_anchored",
    "B_only_anchored",
    "ambiguous_multi_anchor",
    "neither_anchored",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True, help="Normalized JCVI syntenic-pair TSV")
    parser.add_argument("--reference-label", default="gtr")
    parser.add_argument(
        "--group",
        nargs="+",
        required=True,
        metavar="DATASET,SPECIES,A_LABEL,B_LABEL",
        help="One or more A/B groups; PZI and PSP are compared within each dataset.",
    )
    parser.add_argument("--out-states", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-joint", required=True)
    parser.add_argument("--out-candidates", required=True)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_group(text):
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 4 or not all(parts):
        raise ValueError("--group requires DATASET,SPECIES,A_LABEL,B_LABEL: " + text)
    return {"dataset": parts[0], "species": parts[1], "a_label": parts[2], "b_label": parts[3]}


def anchor_state(a_count, b_count):
    if a_count > 1 or b_count > 1:
        return "ambiguous_multi_anchor"
    if a_count == 1 and b_count == 1:
        return "both_anchored"
    if a_count == 1:
        return "A_only_anchored"
    if b_count == 1:
        return "B_only_anchored"
    return "neither_anchored"


def main():
    args = parse_args()
    groups = [parse_group(text) for text in args.group]
    group_by_label = {}
    for group in groups:
        for label, side in ((group["a_label"], "A"), (group["b_label"], "B")):
            if label in group_by_label:
                raise SystemExit("ERROR: query label appears in more than one group: " + label)
            group_by_label[label] = (group, side)

    anchors = defaultdict(lambda: defaultdict(set))
    observed_reference_genes = set()
    for row in read_tsv(args.pairs):
        if row.get("ref_label") != args.reference_label:
            continue
        query_label = row.get("query_label", "")
        if query_label not in group_by_label:
            continue
        ref_gene = row.get("ref_gene", "")
        query_gene = row.get("query_gene", "")
        if not ref_gene or not query_gene:
            continue
        anchors[query_label][ref_gene].add(query_gene)
        observed_reference_genes.add(ref_gene)

    if not observed_reference_genes:
        raise SystemExit(
            "ERROR: no %s anchors found for the requested groups in %s"
            % (args.reference_label, args.pairs)
        )

    state_rows = []
    by_group_and_gene = {}
    for group in groups:
        for ref_gene in sorted(observed_reference_genes):
            a_genes = sorted(anchors[group["a_label"]].get(ref_gene, set()))
            b_genes = sorted(anchors[group["b_label"]].get(ref_gene, set()))
            state = anchor_state(len(a_genes), len(b_genes))
            row = {
                "dataset": group["dataset"],
                "species": group["species"],
                "reference_label": args.reference_label,
                "reference_gene": ref_gene,
                "a_label": group["a_label"],
                "a_anchor_count": len(a_genes),
                "a_anchor_genes": ",".join(a_genes),
                "b_label": group["b_label"],
                "b_anchor_count": len(b_genes),
                "b_anchor_genes": ",".join(b_genes),
                "anchor_state": state,
            }
            state_rows.append(row)
            by_group_and_gene[(group["dataset"], group["species"], ref_gene)] = row

    summary_rows = []
    for group in groups:
        states = [
            by_group_and_gene[(group["dataset"], group["species"], ref_gene)]["anchor_state"]
            for ref_gene in observed_reference_genes
        ]
        counts = Counter(states)
        total = len(states)
        for state in STATE_ORDER:
            count = counts.get(state, 0)
            summary_rows.append(
                {
                    "dataset": group["dataset"],
                    "species": group["species"],
                    "reference_label": args.reference_label,
                    "anchor_state": state,
                    "reference_loci": count,
                    "fraction": "%.8f" % (float(count) / total if total else 0.0),
                    "total_reference_loci": total,
                }
            )

    joint_rows = []
    candidate_rows = []
    datasets = sorted(set(group["dataset"] for group in groups))
    for dataset in datasets:
        pzi = {key[2]: row for key, row in by_group_and_gene.items() if key[0] == dataset and key[1] == "PZI"}
        psp = {key[2]: row for key, row in by_group_and_gene.items() if key[0] == dataset and key[1] == "PSP"}
        if not pzi or not psp:
            continue
        for ref_gene in sorted(observed_reference_genes):
            pzi_row = pzi[ref_gene]
            psp_row = psp[ref_gene]
            left = pzi_row["anchor_state"]
            right = psp_row["anchor_state"]
            if left == "both_anchored" and right == "both_anchored":
                joint_state = "both_anchored_in_both"
            elif left == "A_only_anchored" and right == "A_only_anchored":
                joint_state = "same_A_only_candidate"
            elif left == "B_only_anchored" and right == "B_only_anchored":
                joint_state = "same_B_only_candidate"
            elif {left, right} == {"A_only_anchored", "B_only_anchored"}:
                joint_state = "opposite_single_anchor"
            elif left == right:
                joint_state = "same_ambiguous_state"
            else:
                joint_state = "discordant_or_uninformative"
            joint = {
                "dataset": dataset,
                "reference_label": args.reference_label,
                "reference_gene": ref_gene,
                "pzi_anchor_state": left,
                "psp_anchor_state": right,
                "joint_anchor_state": joint_state,
            }
            joint_rows.append(joint)
            if joint_state in ("same_A_only_candidate", "same_B_only_candidate"):
                candidate = dict(joint)
                candidate.update(
                    {
                        "pzi_A_anchor": pzi_row["a_anchor_genes"],
                        "pzi_B_anchor": pzi_row["b_anchor_genes"],
                        "psp_A_anchor": psp_row["a_anchor_genes"],
                        "psp_B_anchor": psp_row["b_anchor_genes"],
                        "validation_status": "requires_local_synteny_and_sequence_validation",
                    }
                )
                candidate_rows.append(candidate)

    write_tsv(
        args.out_states,
        [
            "dataset", "species", "reference_label", "reference_gene", "a_label",
            "a_anchor_count", "a_anchor_genes", "b_label", "b_anchor_count",
            "b_anchor_genes", "anchor_state",
        ],
        state_rows,
    )
    write_tsv(
        args.out_summary,
        ["dataset", "species", "reference_label", "anchor_state", "reference_loci", "fraction", "total_reference_loci"],
        summary_rows,
    )
    write_tsv(
        args.out_joint,
        ["dataset", "reference_label", "reference_gene", "pzi_anchor_state", "psp_anchor_state", "joint_anchor_state"],
        joint_rows,
    )
    write_tsv(
        args.out_candidates,
        [
            "dataset", "reference_label", "reference_gene", "pzi_anchor_state", "psp_anchor_state",
            "joint_anchor_state", "pzi_A_anchor", "pzi_B_anchor", "psp_A_anchor", "psp_B_anchor",
            "validation_status",
        ],
        candidate_rows,
    )
    print("[OK] Wrote %d anchor states and %d candidate loci." % (len(state_rows), len(candidate_rows)))


if __name__ == "__main__":
    main()
