#!/usr/bin/env python3
"""Concatenate per-locus FASTA alignments and write an IQ-TREE partition file."""

import argparse
from pathlib import Path


def read_fasta(path):
    order = []
    records = {}
    current = None
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].strip().split()[0]
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = "".join(chunks)
    return order, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--pattern", default="*.fa")
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-partitions", required=True)
    args = parser.parse_args()

    concat_order = None
    concat = {}
    partitions = []
    pos = 1
    for path in sorted(Path(args.alignment_dir).glob(args.pattern)):
        order, records = read_fasta(path)
        if concat_order is None:
            concat_order = order
            concat = {taxon: [] for taxon in order}
        if order != concat_order:
            raise SystemExit("ERROR: taxon order mismatch in %s" % path)
        lengths = {len(seq) for seq in records.values()}
        if len(lengths) != 1:
            raise SystemExit("ERROR: unequal sequence lengths in %s" % path)
        length = lengths.pop()
        for taxon in concat_order:
            concat[taxon].append(records[taxon])
        start = pos
        end = pos + length - 1
        partitions.append((path.stem, start, end))
        pos = end + 1

    with Path(args.out_fasta).open("w", encoding="utf-8") as out:
        for taxon in concat_order or []:
            seq = "".join(concat[taxon])
            out.write(">%s\n" % taxon)
            for idx in range(0, len(seq), 80):
                out.write(seq[idx:idx + 80] + "\n")
    with Path(args.out_partitions).open("w", encoding="utf-8") as out:
        for name, start, end in partitions:
            out.write("DNA, %s = %d-%d\n" % (name, start, end))
    print("loci=%d sites=%d" % (len(partitions), pos - 1))


if __name__ == "__main__":
    main()
