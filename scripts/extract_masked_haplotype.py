#!/usr/bin/env python3
import argparse
from pathlib import Path


def read_map(path, copy_name):
    if copy_name not in {"copy1", "copy2"}:
        raise SystemExit("ERROR: copy-name must be copy1 or copy2")
    source_col = "{}_source".format(copy_name)
    mapping = {}
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            normalized_i = header.index("normalized_chr")
            source_i = header.index(source_col)
        except ValueError:
            # The current provenance package stores the corrected pair table
            # instead of the older copy_fastas.map.tsv.  In that table the
            # normalized IDs are copy1_chr/copy2_chr and the source IDs are
            # original_copy1_chr/original_copy2_chr.
            try:
                normalized_i = header.index("{}_chr".format(copy_name))
                source_i = header.index("original_{}_chr".format(copy_name))
            except ValueError as exc:
                raise SystemExit(
                    "ERROR: map file must contain normalized_chr/{0}_source "
                    "or {0}_chr/original_{0}_chr".format(copy_name)
                ) from exc
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            normalized = fields[normalized_i]
            source = fields[source_i]
            mapping[source] = normalized
    if not mapping:
        raise SystemExit("ERROR: no chromosome mappings were read from {}".format(path))
    return mapping


def write_wrapped(out, seq, width=80):
    for i in range(0, len(seq), width):
        out.write(seq[i:i + width] + "\n")


def extract_fasta(input_fasta, mapping, output_fasta):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    found = set()
    current_name = None
    current_seq = []

    def flush(out):
        nonlocal current_name, current_seq
        if current_name is None:
            return
        if current_name in mapping:
            out.write(">{}\n".format(mapping[current_name]))
            write_wrapped(out, "".join(current_seq))
            found.add(current_name)
        current_name = None
        current_seq = []

    with Path(input_fasta).open(encoding="utf-8", errors="replace") as handle, \
            Path(output_fasta).open("w", encoding="utf-8") as out:
        for line in handle:
            if line.startswith(">"):
                flush(out)
                current_name = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        flush(out)

    missing = sorted(set(mapping) - found)
    if missing:
        raise SystemExit("ERROR: missing source chromosomes in masked FASTA: {}".format(",".join(missing)))
    if not output_fasta.is_file() or output_fasta.stat().st_size == 0:
        raise SystemExit("ERROR: empty output FASTA: {}".format(output_fasta))


def parse_args():
    parser = argparse.ArgumentParser(description="Extract a renamed masked copy1/copy2 FASTA from 2n=32 EDTA masked genome.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--map-tsv", required=True)
    parser.add_argument("--copy-name", choices=["copy1", "copy2"], required=True)
    parser.add_argument("--output-fasta", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    mapping = read_map(args.map_tsv, args.copy_name)
    extract_fasta(Path(args.input_fasta), mapping, Path(args.output_fasta))


if __name__ == "__main__":
    main()
