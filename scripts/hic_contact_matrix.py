#!/usr/bin/env python3
import argparse
import collections
import csv
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description="Build a chromosome-pair Hi-C contact matrix from a BAM.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--contact-tsv", required=True)
    parser.add_argument("--matrix-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--heatmap-png", required=True)
    parser.add_argument("--min-mapq", type=int, default=10)
    return parser.parse_args()


def main():
    args = parse_args()

    chroms = []
    with open(args.chrom_sizes) as handle:
        for line in handle:
            if line.strip():
                chroms.append(line.split()[0])
    chrom_set = set(chroms)
    chrom_rank = {chrom: idx for idx, chrom in enumerate(chroms)}

    counts = collections.Counter()
    total_pairs = 0
    used_pairs = 0

    proc = subprocess.Popen(
        ["samtools", "view", "-f", "1", "-F", "3340", args.bam],
        stdout=subprocess.PIPE,
        text=True,
    )

    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        flag = int(fields[1])
        rname = fields[2]
        mapq = int(fields[4])
        mate_rname = fields[6]
        if rname not in chrom_set or mapq < args.min_mapq:
            continue
        if mate_rname == "=":
            mate_rname = rname
        if mate_rname not in chrom_set:
            continue
        if not (flag & 0x40):
            continue
        total_pairs += 1
        a, b = sorted((rname, mate_rname), key=lambda chrom: chrom_rank[chrom])
        counts[(a, b)] += 1
        used_pairs += 1

    return_code = proc.wait()
    if return_code != 0:
        raise SystemExit(f"samtools view failed with exit code {return_code}")

    with open(args.contact_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr_a", "chr_b", "contacts"])
        for i, a in enumerate(chroms):
            for b in chroms[i:]:
                writer.writerow([a, b, counts.get((a, b), 0)])

    with open(args.matrix_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr"] + chroms)
        for a in chroms:
            row = [a]
            for b in chroms:
                key = tuple(sorted((a, b), key=lambda chrom: chrom_rank[chrom]))
                row.append(counts.get(key, 0))
            writer.writerow(row)

    with open(args.summary_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr", "intra_contacts", "inter_contacts", "total_contacts", "cis_percent"])
        for chrom in chroms:
            intra = counts.get((chrom, chrom), 0)
            inter = 0
            for other in chroms:
                if other == chrom:
                    continue
                key = tuple(sorted((chrom, other), key=lambda item: chrom_rank[item]))
                inter += counts.get(key, 0)
            total = intra + inter
            cis_percent = (intra / total * 100) if total else 0
            writer.writerow([chrom, intra, inter, total, f"{cis_percent:.4f}"])

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        matrix = np.zeros((len(chroms), len(chroms)), dtype=float)
        for i, a in enumerate(chroms):
            for j, b in enumerate(chroms):
                key = tuple(sorted((a, b), key=lambda chrom: chrom_rank[chrom]))
                matrix[i, j] = counts.get(key, 0)
        matrix = np.log10(matrix + 1)

        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(matrix, cmap="magma", interpolation="nearest")
        ax.set_xticks(range(len(chroms)))
        ax.set_yticks(range(len(chroms)))
        ax.set_xticklabels(chroms, rotation=90, fontsize=6)
        ax.set_yticklabels(chroms, fontsize=6)
        ax.set_title("Final assembly Hi-C chromosome contact map, log10(count+1)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(args.heatmap_png, dpi=300)
    except Exception as exc:
        with open(args.heatmap_png + ".not_created.txt", "w") as handle:
            handle.write(f"Matplotlib heatmap was not created: {exc}\n")

    print(f"Total read pairs passing filters: {total_pairs}")
    print(f"Paired contacts counted: {used_pairs}")


if __name__ == "__main__":
    main()
