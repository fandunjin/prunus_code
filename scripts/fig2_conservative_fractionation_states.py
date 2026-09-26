#!/usr/bin/env python3
"""Call conservative A/B retention and candidate-fractionation states.

GTR loci are retained only when both PCAM and PPE have exactly one syntenic
anchor.  This is a conservative external-anchor filter, not a verified gene
loss call.  PZI, PSP and PPA are then evaluated as independent A/B groups.
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
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--reference-label", default="gtr")
    parser.add_argument(
        "--polyploid-group", nargs="+", required=True,
        metavar="DATASET,SPECIES,A_LABEL,B_LABEL",
    )
    parser.add_argument(
        "--diploid-control", nargs="+", required=True,
        metavar="SPECIES,LABEL",
        help="Exactly two external controls, normally PCAM and PPE.",
    )
    parser.add_argument("--out-external-status", required=True)
    parser.add_argument("--out-conservative-loci", required=True)
    parser.add_argument("--out-states", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-joint", required=True)
    parser.add_argument("--out-candidates", required=True)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_polyploid(text):
    values = [value.strip() for value in text.split(",")]
    if len(values) != 4 or not all(values):
        raise ValueError("--polyploid-group requires DATASET,SPECIES,A_LABEL,B_LABEL: " + text)
    return {"dataset": values[0], "species": values[1], "a_label": values[2], "b_label": values[3]}


def parse_control(text):
    values = [value.strip() for value in text.split(",")]
    if len(values) != 2 or not all(values):
        raise ValueError("--diploid-control requires SPECIES,LABEL: " + text)
    return {"species": values[0], "label": values[1]}


def state_from_counts(a_count, b_count):
    if a_count > 1 or b_count > 1:
        return "ambiguous_multi_anchor"
    if a_count == 1 and b_count == 1:
        return "both_anchored"
    if a_count == 1:
        return "A_only_anchored"
    if b_count == 1:
        return "B_only_anchored"
    return "neither_anchored"


def joint_state(left, right):
    if left == "both_anchored" and right == "both_anchored":
        return "both_anchored_in_both"
    if left == "A_only_anchored" and right == "A_only_anchored":
        return "same_A_only_candidate"
    if left == "B_only_anchored" and right == "B_only_anchored":
        return "same_B_only_candidate"
    if {left, right} == {"A_only_anchored", "B_only_anchored"}:
        return "opposite_single_anchor"
    if left == right:
        return "same_ambiguous_state"
    return "discordant_or_uninformative"


def main():
    args = parse_args()
    groups = [parse_polyploid(text) for text in args.polyploid_group]
    controls = [parse_control(text) for text in args.diploid_control]
    if len(controls) != 2:
        raise SystemExit("ERROR: exactly two --diploid-control values are required")
    if {control["species"] for control in controls} != {"PCAM", "PPE"}:
        raise SystemExit("ERROR: controls must be named PCAM and PPE")

    group_labels = set()
    for group in groups:
        group_labels.update([group["a_label"], group["b_label"]])
    control_labels = {control["label"] for control in controls}
    requested_labels = group_labels | control_labels

    anchors = defaultdict(lambda: defaultdict(set))
    reference_meta = {}
    observed_reference_genes = set()
    for row in read_tsv(args.pairs):
        if row.get("ref_label") != args.reference_label:
            continue
        label = row.get("query_label", "")
        if label not in requested_labels:
            continue
        ref_gene = row.get("ref_gene", "")
        query_gene = row.get("query_gene", "")
        if not ref_gene or not query_gene:
            continue
        anchors[label][ref_gene].add(query_gene)
        observed_reference_genes.add(ref_gene)
        reference_meta.setdefault(
            ref_gene,
            {
                "reference_seqid": row.get("ref_seqid", ""),
                "reference_start0": row.get("ref_start0", ""),
                "reference_end": row.get("ref_end", ""),
            },
        )
    if not observed_reference_genes:
        raise SystemExit("ERROR: no requested anchors found")

    control_by_species = {control["species"]: control["label"] for control in controls}
    external_rows = []
    conservative_loci = []
    for ref_gene in sorted(observed_reference_genes):
        pcam_genes = sorted(anchors[control_by_species["PCAM"]].get(ref_gene, set()))
        ppe_genes = sorted(anchors[control_by_species["PPE"]].get(ref_gene, set()))
        is_conservative = len(pcam_genes) == 1 and len(ppe_genes) == 1
        status = "conserved_single_copy_external_anchor" if is_conservative else "not_conservative_external_anchor"
        row = {
            "reference_label": args.reference_label,
            "reference_gene": ref_gene,
            "reference_seqid": reference_meta[ref_gene]["reference_seqid"],
            "reference_start0": reference_meta[ref_gene]["reference_start0"],
            "reference_end": reference_meta[ref_gene]["reference_end"],
            "pcam_label": control_by_species["PCAM"],
            "pcam_anchor_count": len(pcam_genes),
            "pcam_anchor_genes": ",".join(pcam_genes),
            "ppe_label": control_by_species["PPE"],
            "ppe_anchor_count": len(ppe_genes),
            "ppe_anchor_genes": ",".join(ppe_genes),
            "external_anchor_status": status,
        }
        external_rows.append(row)
        if is_conservative:
            conservative_loci.append(row)
    if not conservative_loci:
        raise SystemExit("ERROR: no conservative GTR-PCAM-PPE loci remain")

    state_rows = []
    state_by_group = {}
    for group in groups:
        for locus in conservative_loci:
            ref_gene = locus["reference_gene"]
            a_genes = sorted(anchors[group["a_label"]].get(ref_gene, set()))
            b_genes = sorted(anchors[group["b_label"]].get(ref_gene, set()))
            state = state_from_counts(len(a_genes), len(b_genes))
            row = {
                "dataset": group["dataset"],
                "species": group["species"],
                "reference_label": args.reference_label,
                "reference_gene": ref_gene,
                "reference_seqid": locus["reference_seqid"],
                "reference_start0": locus["reference_start0"],
                "reference_end": locus["reference_end"],
                "a_label": group["a_label"],
                "a_anchor_count": len(a_genes),
                "a_anchor_genes": ",".join(a_genes),
                "b_label": group["b_label"],
                "b_anchor_count": len(b_genes),
                "b_anchor_genes": ",".join(b_genes),
                "anchor_state": state,
                "interpretation": "candidate_retention_state_not_verified_fractionation",
            }
            state_rows.append(row)
            state_by_group[(group["dataset"], group["species"], ref_gene)] = row

    summary_rows = []
    for group in groups:
        rows = [
            state_by_group[(group["dataset"], group["species"], locus["reference_gene"])]
            for locus in conservative_loci
        ]
        counts = Counter(row["anchor_state"] for row in rows)
        for state in STATE_ORDER:
            count = counts.get(state, 0)
            summary_rows.append(
                {
                    "dataset": group["dataset"],
                    "species": group["species"],
                    "reference_label": args.reference_label,
                    "anchor_state": state,
                    "reference_loci": count,
                    "fraction": "%.8f" % (float(count) / len(rows)),
                    "total_reference_loci": len(rows),
                    "anchor_universe": "GTR loci with unique PCAM and PPE anchors",
                }
            )

    joint_rows = []
    candidate_rows = []
    for dataset in sorted({group["dataset"] for group in groups}):
        pzi = {row["reference_gene"]: row for row in state_rows if row["dataset"] == dataset and row["species"] == "PZI"}
        psp = {row["reference_gene"]: row for row in state_rows if row["dataset"] == dataset and row["species"] == "PSP"}
        ppa = {row["reference_gene"]: row for row in state_rows if row["dataset"] == dataset and row["species"] == "PPA"}
        if not pzi or not psp:
            continue
        for locus in conservative_loci:
            ref_gene = locus["reference_gene"]
            left = pzi[ref_gene]
            right = psp[ref_gene]
            ppa_row = ppa.get(ref_gene)
            category = joint_state(left["anchor_state"], right["anchor_state"])
            joint = {
                "dataset": dataset,
                "reference_label": args.reference_label,
                "reference_gene": ref_gene,
                "reference_seqid": locus["reference_seqid"],
                "reference_start0": locus["reference_start0"],
                "reference_end": locus["reference_end"],
                "pzi_anchor_state": left["anchor_state"],
                "psp_anchor_state": right["anchor_state"],
                "ppa_anchor_state": ppa_row["anchor_state"] if ppa_row else "not_available",
                "joint_anchor_state": category,
            }
            joint_rows.append(joint)
            if category in {"same_A_only_candidate", "same_B_only_candidate"}:
                candidate = dict(joint)
                candidate.update(
                    {
                        "pzi_A_anchor": left["a_anchor_genes"],
                        "pzi_B_anchor": left["b_anchor_genes"],
                        "psp_A_anchor": right["a_anchor_genes"],
                        "psp_B_anchor": right["b_anchor_genes"],
                        "ppa_A_anchor": ppa_row["a_anchor_genes"] if ppa_row else "",
                        "ppa_B_anchor": ppa_row["b_anchor_genes"] if ppa_row else "",
                        "validation_status": "requires_flanking_synteny_sequence_search_gap_and_read_depth_checks",
                    }
                )
                candidate_rows.append(candidate)

    write_tsv(
        args.out_external_status,
        [
            "reference_label", "reference_gene", "reference_seqid", "reference_start0", "reference_end",
            "pcam_label", "pcam_anchor_count", "pcam_anchor_genes", "ppe_label", "ppe_anchor_count",
            "ppe_anchor_genes", "external_anchor_status",
        ],
        external_rows,
    )
    write_tsv(
        args.out_conservative_loci,
        [
            "reference_label", "reference_gene", "reference_seqid", "reference_start0", "reference_end",
            "pcam_label", "pcam_anchor_count", "pcam_anchor_genes", "ppe_label", "ppe_anchor_count",
            "ppe_anchor_genes", "external_anchor_status",
        ],
        conservative_loci,
    )
    write_tsv(
        args.out_states,
        [
            "dataset", "species", "reference_label", "reference_gene", "reference_seqid", "reference_start0",
            "reference_end", "a_label", "a_anchor_count", "a_anchor_genes", "b_label", "b_anchor_count",
            "b_anchor_genes", "anchor_state", "interpretation",
        ],
        state_rows,
    )
    write_tsv(
        args.out_summary,
        [
            "dataset", "species", "reference_label", "anchor_state", "reference_loci", "fraction",
            "total_reference_loci", "anchor_universe",
        ],
        summary_rows,
    )
    write_tsv(
        args.out_joint,
        [
            "dataset", "reference_label", "reference_gene", "reference_seqid", "reference_start0", "reference_end",
            "pzi_anchor_state", "psp_anchor_state", "ppa_anchor_state", "joint_anchor_state",
        ],
        joint_rows,
    )
    write_tsv(
        args.out_candidates,
        [
            "dataset", "reference_label", "reference_gene", "reference_seqid", "reference_start0", "reference_end",
            "pzi_anchor_state", "psp_anchor_state", "ppa_anchor_state", "joint_anchor_state", "pzi_A_anchor",
            "pzi_B_anchor", "psp_A_anchor", "psp_B_anchor", "ppa_A_anchor", "ppa_B_anchor", "validation_status",
        ],
        candidate_rows,
    )
    print(
        "[OK] Conservative loci=%d; states=%d; PZI/PSP candidate loci=%d."
        % (len(conservative_loci), len(state_rows), len(candidate_rows))
    )


if __name__ == "__main__":
    main()
