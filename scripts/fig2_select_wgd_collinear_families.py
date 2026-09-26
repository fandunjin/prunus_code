#!/usr/bin/env python3
"""Select conservative WGD-family candidates from collinear A-B gene pairs.

This implements the filtering stage used before the Fig. 3b-style gene-family
tree analysis.  It deliberately does not use all OrthoFinder duplication calls:
only orthogroups with one collinear A/B pair in every focal tetraploid and the
specified outgroup configuration are retained.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_id_map(path):
    result = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            original = row["original_id"]
            stable = row["stable_id"]
            # Some source annotations use a compound FASTA ID such as
            # gene-ID|rna-ID, whereas JCVI anchors retain only rna-ID.
            # Both components are aliases of the same frozen stable ID.
            for alias in {original, *[part for part in original.split("|") if part]}:
                if alias in result and result[alias] != stable:
                    raise ValueError("Non-unique original ID in {}: {}".format(path, alias))
                result[alias] = stable
    return result


def parse_key_path(items, argument):
    result = {}
    for item in items:
        if "=" not in item:
            raise ValueError("{} must be LABEL=PATH: {}".format(argument, item))
        label, path = item.split("=", 1)
        if not label or not path:
            raise ValueError("{} must be LABEL=PATH: {}".format(argument, item))
        if label in result:
            raise ValueError("Duplicate {} label: {}".format(argument, label))
        result[label] = Path(path)
    return result


def parse_pair_path(items):
    result = []
    for item in items:
        fields = item.split(",", 2)
        if len(fields) != 3 or not all(fields):
            raise ValueError("--collinear-pair must be LABEL_A,LABEL_B,PATH: {}".format(item))
        result.append((fields[0], fields[1], Path(fields[2])))
    return result


def read_anchor_pairs(path, first_map, second_map):
    pairs = set()
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < 2:
                raise ValueError("Malformed anchor line {} in {}".format(line_number, path))
            left = first_map.get(fields[0])
            right = second_map.get(fields[1])
            if left is None or right is None:
                raise ValueError(
                    "Anchor ID absent from frozen id map at {}:{} ({}, {})".format(
                        path, line_number, fields[0], fields[1]
                    )
                )
            pairs.add((left, right))
    if not pairs:
        raise ValueError("No gene pairs read from {}".format(path))
    return pairs


def read_fasta(path):
    sequences = {}
    current = None
    chunks = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                if current is None:
                    raise ValueError("Sequence before header in {}".format(path))
                chunks.append(line)
    if current is not None:
        sequences[current] = "".join(chunks)
    if not sequences:
        raise ValueError("No sequences in {}".format(path))
    return sequences


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def split_members(value):
    return [] if not value else [member.strip() for member in value.split(",") if member.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--id-map", action="append", default=[], required=True, help="LABEL=ID_MAP.tsv")
    parser.add_argument("--pep", action="append", default=[], required=True, help="LABEL=FROZEN_PEP.fa")
    parser.add_argument(
        "--collinear-pair", action="append", default=[], required=True,
        help="LABEL_A,LABEL_B,ANCHORS.new; provide one for each focal tetraploid",
    )
    parser.add_argument("--outgroup", default="gtr")
    parser.add_argument("--second-outgroup", default="rch")
    parser.add_argument("--second-outgroup-mode", choices=["at_least_one", "exactly_one"], default="exactly_one")
    parser.add_argument(
        "--required-single-copy", action="append", default=[],
        help="Additional taxa that must each have exactly one orthogroup member; repeatable.",
    )
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    map_paths = parse_key_path(args.id_map, "--id-map")
    pep_paths = parse_key_path(args.pep, "--pep")
    pair_specs = parse_pair_path(args.collinear_pair)
    required_single_copy = list(dict.fromkeys(args.required_single_copy))
    if args.outgroup in required_single_copy or args.second_outgroup in required_single_copy:
        raise SystemExit("ERROR: outgroups must not be repeated in --required-single-copy")
    required_labels = {args.outgroup, args.second_outgroup, *required_single_copy}
    for first_label, second_label, _ in pair_specs:
        required_labels.update([first_label, second_label])
    missing_maps = sorted(required_labels - set(map_paths))
    missing_pep = sorted(required_labels - set(pep_paths))
    if missing_maps or missing_pep:
        raise SystemExit("ERROR: missing labels; id-map={}, pep={}".format(missing_maps, missing_pep))

    id_maps = {label: read_id_map(path) for label, path in map_paths.items()}
    anchor_pairs = {}
    for first_label, second_label, path in pair_specs:
        anchor_pairs[(first_label, second_label)] = read_anchor_pairs(
            path, id_maps[first_label], id_maps[second_label]
        )

    selected_rows = []
    all_rows = []
    selected_members = {}
    rejection_counts = Counter()
    with Path(args.orthogroups).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Orthogroup" not in reader.fieldnames:
            raise ValueError("Missing Orthogroup column: {}".format(args.orthogroups))
        missing_columns = sorted(required_labels - set(reader.fieldnames))
        if missing_columns:
            raise ValueError("Orthogroups lacks required labels: {}".format(",".join(missing_columns)))
        for row in reader:
            og = row["Orthogroup"]
            members = {label: split_members(row.get(label, "")) for label in required_labels}
            reasons = []
            pair_values = {}
            for first_label, second_label, _ in pair_specs:
                first_members = members[first_label]
                second_members = members[second_label]
                key = first_label + "__" + second_label
                if len(first_members) != 1 or len(second_members) != 1:
                    reasons.append(key + "_not_one_copy_each")
                    pair_values[key] = ""
                    continue
                pair = (first_members[0], second_members[0])
                if pair not in anchor_pairs[(first_label, second_label)]:
                    reasons.append(key + "_not_collinear")
                    pair_values[key] = ""
                else:
                    pair_values[key] = "{}|{}".format(pair[0], pair[1])
            if len(members[args.outgroup]) != 1:
                reasons.append(args.outgroup + "_not_one_copy")
            if args.second_outgroup_mode == "exactly_one" and len(members[args.second_outgroup]) != 1:
                reasons.append(args.second_outgroup + "_not_one_copy")
            elif args.second_outgroup_mode == "at_least_one" and len(members[args.second_outgroup]) < 1:
                reasons.append(args.second_outgroup + "_absent")
            for label in required_single_copy:
                if len(members[label]) != 1:
                    reasons.append(label + "_not_one_copy")
            eligible = not reasons
            if eligible:
                selected_rows.append({
                    "orthogroup": og,
                    **pair_values,
                    args.outgroup + "_copy": members[args.outgroup][0],
                    args.second_outgroup + "_copies": ",".join(members[args.second_outgroup]),
                    **{label + "_copy": members[label][0] for label in required_single_copy},
                })
                selected_members[og] = members
            else:
                for reason in reasons:
                    rejection_counts[reason] += 1
            all_rows.append({"orthogroup": og, "eligible": "yes" if eligible else "no", "rejection_reasons": ";".join(reasons), **pair_values})

    fields = ["orthogroup"] + [first + "__" + second for first, second, _ in pair_specs] + [
        args.outgroup + "_copy", args.second_outgroup + "_copies"
    ] + [label + "_copy" for label in required_single_copy]
    write_tsv(outdir / "selected_collinear_wgd_families.tsv", fields, selected_rows)
    write_tsv(outdir / "all_orthogroup_screen.tsv", ["orthogroup", "eligible", "rejection_reasons"] + fields[1:-2], all_rows)
    write_tsv(
        outdir / "screening_summary.tsv",
        ["metric", "count"],
        [{"metric": "orthogroups_screened", "count": len(all_rows)}, {"metric": "families_retained", "count": len(selected_rows)}]
        + [{"metric": "rejected_" + key, "count": value} for key, value in sorted(rejection_counts.items())],
    )

    fasta_by_label = {label: read_fasta(pep_paths[label]) for label in required_labels}
    fasta_dir = outdir / "family_fastas"
    fasta_dir.mkdir()
    member_rows = []
    for og, members in selected_members.items():
        chosen = []
        for first_label, second_label, _ in pair_specs:
            chosen.extend([(first_label, members[first_label][0]), (second_label, members[second_label][0])])
        chosen.append((args.outgroup, members[args.outgroup][0]))
        chosen.extend((args.second_outgroup, member) for member in members[args.second_outgroup])
        chosen.extend((label, members[label][0]) for label in required_single_copy)
        # A designated outgroup can also occur in a cross-species collinear
        # pair (for example PCAM-GTR).  It must occur only once in each tree.
        unique_chosen = []
        seen = set()
        for item in chosen:
            if item not in seen:
                unique_chosen.append(item)
                seen.add(item)
        with (fasta_dir / (og + ".fa")).open("w", encoding="utf-8", newline="\n") as handle:
            for label, stable_id in unique_chosen:
                sequence = fasta_by_label[label].get(stable_id)
                if sequence is None:
                    raise ValueError("Selected ID absent from {} PEP: {}".format(label, stable_id))
                handle.write(">{0}|{1}\n{2}\n".format(label, stable_id, sequence))
                member_rows.append({"orthogroup": og, "label": label, "stable_id": stable_id})
    write_tsv(outdir / "selected_family_members.tsv", ["orthogroup", "label", "stable_id"], member_rows)
    (outdir / "method_note.txt").write_text(
        "Selection follows the Fig. 3b pre-tree filter: each retained family has one A-B collinear pair in every focal tetraploid, one GTR copy, {} RCH copy/copies, and one copy in each additional required taxon ({}).\n"
        "This is an input-selection stage only. Duplication placement requires subsequent family-tree reconstruction and reconciliation/PUG.\n".format(
            "exactly one" if args.second_outgroup_mode == "exactly_one" else "at least one",
            ",".join(required_single_copy) if required_single_copy else "none",
        ),
        encoding="ascii",
    )
    print("[OK] Screened {} orthogroups; retained {} conservative collinear WGD-family candidates.".format(
        len(all_rows), len(selected_rows)
    ))


if __name__ == "__main__":
    main()
