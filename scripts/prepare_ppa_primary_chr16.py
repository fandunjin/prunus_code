#!/usr/bin/env python3
"""Record the 16 chromosome IDs in the annotated PPA primary genome.

An optional neutral FASTA export is available for inspection, but Figure 2
minimap2 alignment should use the original primary assembly itself. Biological
PPA_A/PPA_B names are applied later, after an explicitly reviewed
PPA-to-PZI correspondence table is available.
"""

import argparse
import csv
import gzip
import re
from pathlib import Path


CHROMOSOME_RE = re.compile(r"\bchromosome\s*[:=]?\s*(\d+)\b", re.IGNORECASE)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def fasta_records(path: Path):
    name = None
    header = ""
    chunks = []
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = raw.rstrip("\r\n")[1:]
                name = header.split(None, 1)[0]
                chunks = []
            else:
                chunks.append(raw.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--out-fasta", type=Path)
    parser.add_argument("--out-map", type=Path, required=True)
    args = parser.parse_args()

    if not args.genome.is_file() or args.genome.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing primary genome: {args.genome}")
    records = []
    seen_numbers = set()
    seen_names = set()
    for source_id, header, sequence in fasta_records(args.genome):
        match = CHROMOSOME_RE.search(header)
        if match is None:
            continue
        number = int(match.group(1))
        if number in seen_numbers:
            raise SystemExit(f"ERROR: duplicate chromosome number {number} in {args.genome}")
        if source_id in seen_names:
            raise SystemExit(f"ERROR: duplicate FASTA record ID {source_id}")
        seen_numbers.add(number)
        seen_names.add(source_id)
        records.append((number, source_id, header, sequence))

    records.sort(key=lambda row: row[0])
    expected = list(range(1, 17))
    observed = [row[0] for row in records]
    if observed != expected:
        raise SystemExit(f"ERROR: expected chromosome labels 1--16, observed {observed}")

    args.out_map.parent.mkdir(parents=True, exist_ok=True)
    out_fasta = None
    if args.out_fasta is not None:
        args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
        out_fasta = args.out_fasta.open("w", encoding="utf-8", newline="\n")
    with args.out_map.open("w", encoding="utf-8", newline="") as out_map:
        writer = csv.writer(out_map, delimiter="\t", lineterminator="\n")
        writer.writerow(["ppa_source_seqid", "ppa_neutral_seqid", "chromosome_number", "length_bp", "source_header"])
        for number, source_id, header, sequence in records:
            neutral_id = f"PPA_chr{number:02d}"
            if out_fasta is not None:
                out_fasta.write(f">{neutral_id} source_seqid={source_id} chromosome={number}\n")
                for start in range(0, len(sequence), 80):
                    out_fasta.write(sequence[start:start + 80] + "\n")
            writer.writerow([source_id, neutral_id, number, len(sequence), header])
    if out_fasta is not None:
        out_fasta.close()
    print(f"[OK] Recorded {len(records)} PPA primary chromosome IDs: {args.out_map}")


if __name__ == "__main__":
    main()
