#!/usr/bin/env python3
"""Prepare the primary-only input manifest for Fig. 2c FraBias analysis.

The script does not run LAST or FraBias.  It freezes the eight tracks and
reports all GTR chromosome anchors for the three selected PZI A/B chromosome
pairs before computationally expensive whole-genome synteny is started.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path


TRACKS = [
    ("PZI_A", "pzi_hapA1", "PZI", "A"),
    ("PZI_B", "pzi_hapB1", "PZI", "B"),
    ("PSP_A", "psp_hapA1", "PSP", "A"),
    ("PSP_B", "psp_hapB1", "PSP", "B"),
    ("PPA_A", "ppa_hapA", "PPA", "A"),
    ("PPA_B", "ppa_hapB", "PPA", "B"),
    ("PCAM", "pcam", "PCAM", "single_copy_control"),
    ("PPE", "ppe", "PPE", "single_copy_control"),
]

PZI_PAIRS = [
    ("C06", "C06_A1", "C06_B1"),
    ("C02", "C02_A1", "C02_B1"),
    ("C03", "C03_A1", "C03_B1"),
]


def chrom_key(value):
    digits = "".join(character for character in value if character.isdigit())
    return int(digits) if digits else 99999


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def count_records(path):
    count = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jcvi-inputs", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    input_root = Path(args.jcvi_inputs)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    track_rows = []
    for track_id, label, species, copy_state in TRACKS:
        source = input_root / label
        bed = source / (label + ".bed")
        pep = source / (label + ".pep")
        seqids = source / (label + ".seqids")
        for path in [bed, pep, seqids]:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError("Missing input: {}".format(path))
        track_rows.append({
            "track_id": track_id,
            "source_label": label,
            "species": species,
            "copy_state": copy_state,
            "bed": str(bed.resolve()),
            "pep": str(pep.resolve()),
            "seqids": str(seqids.resolve()),
            "pep_records": count_records(pep),
            "renamed_gene_prefix": track_id + "__",
            "renamed_chromosome_prefix": track_id + "__",
        })
    write_tsv(
        outdir / "eight_track_manifest.tsv",
        list(track_rows[0].keys()),
        track_rows,
    )

    rows = read_tsv(args.pairs)
    required = {"query_label", "query_seqid", "ref_label", "ref_seqid"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Pair table lacks required columns: {}".format(args.pairs))
    selected = {"pzi_hapA1": {pair[1]: pair[0] for pair in PZI_PAIRS}, "pzi_hapB1": {pair[2]: pair[0] for pair in PZI_PAIRS}}
    counts = Counter()
    for row in rows:
        label = row["query_label"]
        if label not in selected or row["ref_label"] != "gtr":
            continue
        pair_name = selected[label].get(row["query_seqid"])
        if pair_name:
            counts[(pair_name, label, row["query_seqid"], row["ref_seqid"])] += 1

    chromosome_rows = []
    for pair_name, a_chromosome, b_chromosome in PZI_PAIRS:
        for label, chromosome in [("pzi_hapA1", a_chromosome), ("pzi_hapB1", b_chromosome)]:
            values = [
                (ref, count)
                for (this_pair, this_label, this_chromosome, ref), count in counts.items()
                if this_pair == pair_name and this_label == label and this_chromosome == chromosome
            ]
            values.sort(key=lambda item: (-item[1], chrom_key(item[0]), item[0]))
            total = sum(count for _, count in values)
            if not values:
                chromosome_rows.append({
                    "pzi_pair": pair_name,
                    "pzi_label": label,
                    "pzi_haplotype_chromosome": chromosome,
                    "gtr_chromosome": "",
                    "anchor_count": 0,
                    "fraction_of_selected_chromosome_anchors": "0.000000",
                    "rank": 0,
                })
            else:
                for rank, (ref, count) in enumerate(values, start=1):
                    chromosome_rows.append({
                        "pzi_pair": pair_name,
                        "pzi_label": label,
                        "pzi_haplotype_chromosome": chromosome,
                        "gtr_chromosome": ref,
                        "anchor_count": count,
                        "fraction_of_selected_chromosome_anchors": "{:.6f}".format(count / float(total)),
                        "rank": rank,
                    })
    write_tsv(
        outdir / "selected_pzi_chromosome_to_all_gtr_chromosomes.tsv",
        ["pzi_pair", "pzi_label", "pzi_haplotype_chromosome", "gtr_chromosome", "anchor_count", "fraction_of_selected_chromosome_anchors", "rank"],
        chromosome_rows,
    )
    write_tsv(
        outdir / "selected_pzi_homoeologous_pairs.tsv",
        ["pair_label", "pzi_A1_chromosome", "pzi_B1_chromosome"],
        [
            {"pair_label": pair_name, "pzi_A1_chromosome": a_chromosome, "pzi_B1_chromosome": b_chromosome}
            for pair_name, a_chromosome, b_chromosome in PZI_PAIRS
        ],
    )
    (outdir / "method_note.txt").write_text(
        "This primary-only manifest defines the eight Fig. 2c tracks and reports all GTR chromosome anchors for PZI C06/C02/C03 A1-B1 pairs.\n"
        "The subsequent eight-track calculation remains LAST, Quota Align Merge (20 genes), Quota Align (8:1), then FraBias (100 genes).\n"
        "Quota 8:1 is the project-specific depth for PZI_A/PZI_B/PSP_A/PSP_B/PPA_A/PPA_B/PCAM/PPE against GTR.\n"
        "The anchor counts here select and audit display regions; they are not FraBias retention estimates.\n",
        encoding="ascii",
    )
    print("[OK] Wrote primary-only Fig. 2c FraBias manifest: {}".format(outdir))


if __name__ == "__main__":
    main()
