#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_fasta(path):
    name = None
    chunks = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def read_paf_summary(path, chr_count):
    summary = {
        "chr{:02d}".format(i): {"+": 0, "-": 0, "other_chr": 0}
        for i in range(1, chr_count + 1)
    }
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            strand = fields[4]
            tname = fields[5]
            aln_len = int(fields[10])
            if qname not in summary:
                continue
            if tname == qname and strand in {"+", "-"}:
                summary[qname][strand] += aln_len
            else:
                summary[qname]["other_chr"] += aln_len
    return summary


def summarize(args):
    summary = read_paf_summary(args.paf, args.chr_count)
    with Path(args.out_tsv).open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chr",
            "same_chr_plus_aligned_bp",
            "same_chr_minus_aligned_bp",
            "other_chr_aligned_bp",
            "dominant_orientation",
            "same_chr_aligned_bp",
            "same_chr_orientation_ratio",
        ])
        for i in range(1, args.chr_count + 1):
            chrom = "chr{:02d}".format(i)
            plus = summary[chrom]["+"]
            minus = summary[chrom]["-"]
            same = plus + minus
            dominant = "+" if plus >= minus else "-"
            ratio = 0 if same == 0 else max(plus, minus) / same
            writer.writerow([
                chrom,
                plus,
                minus,
                summary[chrom]["other_chr"],
                dominant,
                same,
                "{:.6f}".format(ratio),
            ])


def orient_query(args):
    summary = read_paf_summary(args.paf, args.chr_count)
    orientations = {}
    with Path(args.orientation_tsv).open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_chr",
            "target_chr",
            "orientation",
            "same_chr_plus_aligned_bp",
            "same_chr_minus_aligned_bp",
            "same_chr_aligned_bp",
            "other_chr_aligned_bp",
        ])
        for i in range(1, args.chr_count + 1):
            chrom = "chr{:02d}".format(i)
            plus = summary[chrom]["+"]
            minus = summary[chrom]["-"]
            orientation = "+" if plus >= minus else "-"
            orientations[chrom] = orientation
            writer.writerow([
                chrom,
                chrom,
                orientation,
                plus,
                minus,
                plus + minus,
                summary[chrom]["other_chr"],
            ])

    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    seen = set()
    with Path(args.query_out).open("w", encoding="utf-8") as out:
        for name, seq in read_fasta(args.query_in):
            if name not in orientations:
                continue
            seen.add(name)
            if orientations[name] == "-":
                seq = seq.translate(comp)[::-1]
            out.write(">{}\n".format(name))
            for line in wrap(seq):
                out.write(line + "\n")
    expected = {"chr{:02d}".format(i) for i in range(1, args.chr_count + 1)}
    missing = sorted(expected - seen)
    if missing:
        raise SystemExit("ERROR: missing query chromosomes: {}".format(",".join(missing)))


def parse_args():
    parser = argparse.ArgumentParser(description="Helpers for copy1-vs-copy2 haplotype SyRI workflow.")
    subparsers = parser.add_subparsers(dest="command")

    p1 = subparsers.add_parser("summarize-paf")
    p1.add_argument("--paf", required=True)
    p1.add_argument("--out-tsv", required=True)
    p1.add_argument("--chr-count", type=int, default=16)
    p1.set_defaults(func=summarize)

    p2 = subparsers.add_parser("orient-query")
    p2.add_argument("--paf", required=True)
    p2.add_argument("--query-in", required=True)
    p2.add_argument("--query-out", required=True)
    p2.add_argument("--orientation-tsv", required=True)
    p2.add_argument("--chr-count", type=int, default=16)
    p2.set_defaults(func=orient_query)

    return parser.parse_args()


def main():
    args = parse_args()
    if not hasattr(args, "func"):
        raise SystemExit("ERROR: a subcommand is required.")
    args.func(args)


if __name__ == "__main__":
    main()
