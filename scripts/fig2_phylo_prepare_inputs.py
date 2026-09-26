#!/usr/bin/env python3
"""Freeze and normalize protein FASTA inputs for Figure 2 phylogenomics."""

import argparse
import csv
import hashlib
import os
from pathlib import Path


AA = set("ABCDEFGHIKLMNPQRSTVWXYZJUO")


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path):
    header = None
    chunks = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError("sequence before first FASTA header")
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalize_one(row, destination, map_path, skip_invalid=False):
    source = Path(row["source_pep"])
    if not source.is_file() or source.stat().st_size == 0:
        raise RuntimeError("missing or empty protein FASTA: %s" % source)

    label = row["label"]
    seen = set()
    mapping = []
    record_count = 0
    aa_count = 0
    min_len = None
    max_len = 0
    excluded = []

    # Write to temporary files so a failed freeze cannot leave a partial FASTA
    # that looks like a valid immutable output on the next run.
    destination_tmp = destination.with_name(destination.name + ".tmp")
    map_path_tmp = map_path.with_name(map_path.name + ".tmp")

    with destination_tmp.open("w", encoding="ascii", newline="\n") as output:
        for header, raw_sequence in fasta_records(source):
            original_id = header.split()[0]
            if not original_id:
                raise RuntimeError("empty FASTA identifier in %s" % source)
            if original_id in seen:
                raise RuntimeError("duplicate FASTA identifier %s in %s" % (original_id, source))
            seen.add(original_id)

            sequence = raw_sequence.replace(" ", "").upper()
            sequence = sequence.rstrip("*")
            if not sequence:
                raise RuntimeError("empty protein sequence for %s in %s" % (original_id, source))
            if "*" in sequence:
                positions = [str(index + 1) for index, symbol in enumerate(sequence) if symbol == "*"]
                if not skip_invalid:
                    raise RuntimeError("internal stop codon for %s in %s" % (original_id, source))
                excluded.append({
                    "dataset": row["dataset"],
                    "label": label,
                    "source_pep": str(source),
                    "original_id": original_id,
                    "original_header": header,
                    "reason": "internal_stop",
                    "positions_1based": ",".join(positions),
                })
                continue
            bad = sorted(set(sequence) - AA)
            if bad:
                raise RuntimeError("invalid amino-acid symbols %s for %s in %s" % ("".join(bad), original_id, source))

            record_count += 1
            stable_id = "%s__g%08d" % (label, record_count)
            output.write(">%s\n" % stable_id)
            for start in range(0, len(sequence), 60):
                output.write(sequence[start:start + 60] + "\n")
            mapping.append({
                "stable_id": stable_id,
                "original_id": original_id,
                "original_header": header,
                "protein_length": str(len(sequence)),
            })
            aa_count += len(sequence)
            min_len = len(sequence) if min_len is None else min(min_len, len(sequence))
            max_len = max(max_len, len(sequence))

    if record_count == 0:
        raise RuntimeError("no valid protein records in %s" % source)
    write_tsv(map_path_tmp, mapping, ["stable_id", "original_id", "original_header", "protein_length"])
    destination_tmp.replace(destination)
    map_path_tmp.replace(map_path)
    return {
        "dataset": row["dataset"],
        "label": label,
        "scientific_name": row["scientific_name"],
        "role": row["role"],
        "source_pep": str(source),
        "source_sha256": sha256(source),
        "source_bytes": str(source.stat().st_size),
        "records": str(record_count),
        "amino_acids": str(aa_count),
        "min_length": str(min_len),
        "max_length": str(max_len),
        "normalized_pep": str(destination),
        "id_map": str(map_path),
        "isoform_policy": row["isoform_policy"],
        "status": "PASS",
    }, excluded


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--skip-invalid-proteins",
        action="store_true",
        help="Exclude records containing internal stop codons and record them in excluded_invalid_proteins.tsv",
    )
    args = parser.parse_args()

    manifest = Path(args.manifest).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(manifest)
    if not rows:
        raise SystemExit("ERROR: empty input manifest: %s" % manifest)

    expanded = []
    for row in rows:
        datasets = [row["dataset"]] if row["dataset"] != "both" else ["main", "sensitivity"]
        for dataset in datasets:
            copy = dict(row)
            copy["dataset"] = dataset
            expanded.append(copy)

    keys = set()
    for row in expanded:
        key = (row["dataset"], row["label"])
        if key in keys:
            raise SystemExit("ERROR: duplicate dataset/label in manifest: %s/%s" % key)
        keys.add(key)

    inventory = []
    exclusions = []
    excluded_path = out_dir / "excluded_invalid_proteins.tsv"
    for row in expanded:
        dataset_dir = out_dir / row["dataset"]
        pep_dir = dataset_dir / "pep"
        map_dir = dataset_dir / "id_maps"
        pep_dir.mkdir(parents=True, exist_ok=True)
        map_dir.mkdir(parents=True, exist_ok=True)
        destination = pep_dir / (row["label"] + ".fa")
        map_path = map_dir / (row["label"] + ".id_map.tsv")
        if destination.exists() or map_path.exists():
            raise SystemExit("ERROR: frozen output already exists: %s" % destination)
        result, row_exclusions = normalize_one(
            row,
            destination,
            map_path,
            skip_invalid=args.skip_invalid_proteins,
        )
        inventory.append(result)
        exclusions.extend(row_exclusions)

    fields = [
        "dataset", "label", "scientific_name", "role", "source_pep", "source_sha256",
        "source_bytes", "records", "amino_acids", "min_length", "max_length",
        "normalized_pep", "id_map", "isoform_policy", "status",
    ]
    write_tsv(out_dir / "input_inventory.tsv", inventory, fields)
    write_tsv(excluded_path, exclusions, [
        "dataset", "label", "source_pep", "original_id", "original_header",
        "reason", "positions_1based",
    ])
    with (out_dir / "frozen_input_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(expanded[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(expanded)
    (out_dir / "COMPLETE").write_text("Figure 2 input freeze completed successfully.\n", encoding="ascii")


if __name__ == "__main__":
    main()
