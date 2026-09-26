#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_attrs(text):
    attrs = {}
    for item in text.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(records, path):
    with Path(path).open("w", encoding="utf-8") as out:
        for name, seq in records.items():
            out.write(">{}\n".format(name))
            for idx in range(0, len(seq), 60):
                out.write(seq[idx:idx + 60] + "\n")


def parse_qc_passes(path):
    selected = {}
    order = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            flags = row.get("flags", "")
            if flags != "PASS":
                continue
            gene_id = row["gene_id"]
            mrna_id = row["mrna_id"]
            selected[gene_id] = row
            selected[mrna_id] = row
            order.append((gene_id, mrna_id, row))
    return selected, order


def collect_records(gff3_path, selected_ids):
    lines = []
    feature_counts = Counter()
    seq_counts = Counter()
    kept_gene_ids = set()
    kept_mrna_ids = set()

    with Path(gff3_path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs_text = fields
            attrs = parse_attrs(attrs_text)
            feature_counts[feature] += 1
            keep = False
            if feature == "gene":
                gene_id = attrs.get("ID")
                if gene_id in selected_ids:
                    keep = True
                    kept_gene_ids.add(gene_id)
            elif feature == "mRNA":
                mrna_id = attrs.get("ID")
                parent = attrs.get("Parent", "")
                if mrna_id in selected_ids or parent in selected_ids:
                    keep = True
                    if mrna_id:
                        kept_mrna_ids.add(mrna_id)
            else:
                parent = attrs.get("Parent", "")
                if parent:
                    parents = parent.split(",")
                    if any(p in selected_ids or p in kept_mrna_ids for p in parents):
                        keep = True
            if keep:
                lines.append(line)
                seq_counts[seqid] += 1
    return lines, feature_counts, seq_counts, kept_gene_ids, kept_mrna_ids


def summarize_lengths(records):
    lengths = [len(seq) for seq in records.values()]
    if not lengths:
        return {"count": 0, "min": 0, "median": 0, "mean": 0, "max": 0}
    lengths.sort()
    n = len(lengths)
    median = lengths[n // 2] if n % 2 else (lengths[n // 2 - 1] + lengths[n // 2]) / 2.0
    return {
        "count": n,
        "min": lengths[0],
        "median": median,
        "mean": sum(lengths) / n,
        "max": lengths[-1],
    }


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    selected, ordered = parse_qc_passes(args.qc_tsv)
    if not ordered:
        raise SystemExit("ERROR: no PASS models found in {}".format(args.qc_tsv))

    protein = read_fasta(args.protein_fa)
    cds = read_fasta(args.cds_fa)
    mrna = read_fasta(args.mrna_fa)

    selected_mrna_ids = {row["mrna_id"] for _, _, row in ordered}
    filtered_protein = {k: v for k, v in protein.items() if k in selected_mrna_ids}
    filtered_cds = {k: v for k, v in cds.items() if k in selected_mrna_ids}
    filtered_mrna = {k: v for k, v in mrna.items() if k in selected_mrna_ids}

    kept_lines, feature_counts, seq_counts, kept_gene_ids, kept_mrna_ids = collect_records(
        args.gff3, selected
    )

    filtered_gff3 = out_dir / "{}.high_confidence.gff3".format(args.prefix)
    filtered_protein_fa = out_dir / "{}.high_confidence.protein.fa".format(args.prefix)
    filtered_cds_fa = out_dir / "{}.high_confidence.cds.fa".format(args.prefix)
    filtered_mrna_fa = out_dir / "{}.high_confidence.mrna.fa".format(args.prefix)
    pass_qc_tsv = out_dir / "{}.high_confidence.qc.tsv".format(args.prefix)
    selected_tsv = out_dir / "{}.high_confidence.selected_ids.tsv".format(args.prefix)
    summary_tsv = out_dir / "{}.high_confidence.summary.tsv".format(args.prefix)
    seq_tsv = out_dir / "{}.high_confidence.genes_by_sequence.tsv".format(args.prefix)

    with filtered_gff3.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for line in kept_lines:
            out.write(line if line.endswith("\n") else line + "\n")

    write_fasta(filtered_protein, filtered_protein_fa)
    write_fasta(filtered_cds, filtered_cds_fa)
    write_fasta(filtered_mrna, filtered_mrna_fa)

    with pass_qc_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=ordered[0][2].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for _, _, row in ordered:
            writer.writerow(row)

    with selected_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "mrna_id"])
        for gene_id, mrna_id, _ in ordered:
            writer.writerow([gene_id, mrna_id])

    with seq_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["seq_id", "gene_count"])
        for seqid in sorted(seq_counts):
            writer.writerow([seqid, seq_counts[seqid]])

    protein_stats = summarize_lengths(filtered_protein)
    cds_stats = summarize_lengths(filtered_cds)
    mrna_stats = summarize_lengths(filtered_mrna)

    with summary_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["pass_gene_models", len(ordered)])
        writer.writerow(["input_qc_rows", len(selected_mrna_ids)])
        writer.writerow(["retained_gene_gff3_features", len(kept_lines)])
        writer.writerow(["retained_gene_ids", len(kept_gene_ids)])
        writer.writerow(["retained_mrna_ids", len(kept_mrna_ids)])
        writer.writerow(["retained_protein_records", len(filtered_protein)])
        writer.writerow(["retained_cds_records", len(filtered_cds)])
        writer.writerow(["retained_mrna_records", len(filtered_mrna)])
        for feature, count in sorted(feature_counts.items()):
            writer.writerow(["input_feature_{}".format(feature), count])
        for seqid, count in sorted(seq_counts.items()):
            writer.writerow(["seq_{}_genes".format(seqid), count])
        for name, stats in [("protein", protein_stats), ("cds", cds_stats), ("mrna", mrna_stats)]:
            for key, value in stats.items():
                writer.writerow(["{}_{}".format(name, key), value])

    print("High-confidence set extraction completed")
    print("GFF3:", filtered_gff3)
    print("Protein:", filtered_protein_fa)
    print("CDS:", filtered_cds_fa)
    print("mRNA:", filtered_mrna_fa)
    print("QC:", pass_qc_tsv)
    print("Summary:", summary_tsv)


def parse_args():
    parser = argparse.ArgumentParser(description="Extract PASS high-confidence Helixer gene set.")
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--qc-tsv", required=True)
    parser.add_argument("--protein-fa", required=True)
    parser.add_argument("--cds-fa", required=True)
    parser.add_argument("--mrna-fa", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    main()
