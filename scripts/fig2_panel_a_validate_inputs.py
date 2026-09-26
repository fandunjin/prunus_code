#!/usr/bin/env python3
"""Validate the frozen protein/CDS/BED tracks for strict Figure 2a."""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


TRACKS = (
    "PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B",
    "PCAM", "PPE", "GTR",
)


def read_fasta(path):
    ids = []
    lengths = {}
    current = None
    length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].split(None, 1)[0]
                if not current:
                    raise ValueError("Empty FASTA identifier at {}:{}".format(path, line_number))
                if current in lengths or current in ids:
                    raise ValueError("Duplicate FASTA identifier {} in {}".format(current, path))
                ids.append(current)
                length = 0
            else:
                if current is None:
                    raise ValueError("Sequence before first header in {}".format(path))
                length += len(line)
    if current is not None:
        lengths[current] = length
    if not ids:
        raise ValueError("No FASTA records in {}".format(path))
    return ids, lengths


def read_bed(path):
    ids = []
    seen = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                raise ValueError("Malformed BED at {}:{}".format(path, line_number))
            gene = fields[3]
            if gene in seen:
                raise ValueError("Duplicate BED gene {} in {}".format(gene, path))
            seen.add(gene)
            ids.append(gene)
    if not ids:
        raise ValueError("No BED records in {}".format(path))
    return ids


def write_tsv(path, rows):
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tracks", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    failed = False
    for track in TRACKS:
        pep_path = args.tracks / (track + ".pep")
        cds_path = args.tracks / (track + ".cds")
        bed_path = args.tracks / (track + ".bed")
        for path in (pep_path, cds_path, bed_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit("ERROR: missing input: {}".format(path))

        pep_ids, pep_lengths = read_fasta(pep_path)
        cds_ids, cds_lengths = read_fasta(cds_path)
        bed_ids = read_bed(bed_path)
        pep_set = set(pep_ids)
        cds_set = set(cds_ids)
        bed_set = set(bed_ids)
        pep_cds_match = pep_set == cds_set
        bed_unknown = len(bed_set - pep_set)
        cds_not_triplet = sum(1 for value in cds_lengths.values() if value % 3 != 0)
        empty_pep = sum(1 for value in pep_lengths.values() if value == 0)
        empty_cds = sum(1 for value in cds_lengths.values() if value == 0)
        status = "PASS"
        if not pep_cds_match or bed_unknown or empty_pep or empty_cds:
            status = "FAIL"
            failed = True
        rows.append({
            "track": track,
            "pep_records": len(pep_ids),
            "cds_records": len(cds_ids),
            "bed_records": len(bed_ids),
            "pep_cds_id_match": "yes" if pep_cds_match else "no",
            "bed_ids_absent_from_fasta": bed_unknown,
            "cds_length_not_multiple_of_3": cds_not_triplet,
            "empty_pep_records": empty_pep,
            "empty_cds_records": empty_cds,
            "status": status,
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out, rows)
    print("[OK] Validated {} Figure 2a tracks: {}".format(len(rows), args.out))
    if failed:
        raise SystemExit("ERROR: one or more Figure 2a tracks failed identifier validation")


if __name__ == "__main__":
    main()
