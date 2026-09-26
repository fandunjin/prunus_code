#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def wrap(seq, width=60):
    for idx in range(0, len(seq), width):
        yield seq[idx:idx + width]


def revcomp(seq):
    return seq.translate(COMP)[::-1].upper()


def extract_parts(genome, seqid, intervals, strand):
    seq = genome.get(seqid)
    if seq is None:
        return ""
    ordered = sorted(intervals, key=lambda item: item[0])
    pieces = []
    for start, end in ordered:
        pieces.append(seq[start - 1:end])
    out = "".join(pieces).upper()
    if strand == "-":
        out = revcomp(out)
    return out


def translate(cds):
    aa = []
    upper = cds.upper()
    for idx in range(0, len(upper) - 2, 3):
        codon = upper[idx:idx + 3]
        aa.append(GENETIC_CODE.get(codon, "X"))
    return "".join(aa)


def parse_gff3(path):
    genes = {}
    mrnas = {}
    children = defaultdict(lambda: defaultdict(list))
    feature_counts = Counter()
    duplicate_ids = 0
    ids = set()
    parent_refs = set()
    bad_coords = 0

    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs_text = fields
            feature_counts[feature] += 1
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                bad_coords += 1
                continue
            if start_i < 1 or end_i < start_i:
                bad_coords += 1
            attrs = parse_attrs(attrs_text)
            record = {
                "seqid": seqid,
                "source": source,
                "feature": feature,
                "start": start_i,
                "end": end_i,
                "score": score,
                "strand": strand,
                "phase": phase,
                "attrs": attrs,
                "raw": line.rstrip("\n"),
            }
            rec_id = attrs.get("ID")
            if rec_id:
                if rec_id in ids:
                    duplicate_ids += 1
                ids.add(rec_id)
            parent = attrs.get("Parent")
            if parent:
                for parent_id in parent.split(","):
                    parent_refs.add(parent_id)
                    children[parent_id][feature].append(record)
            if feature == "gene" and rec_id:
                genes[rec_id] = record
            elif feature == "mRNA" and rec_id:
                mrnas[rec_id] = record

    missing_parent_refs = sum(1 for parent_id in parent_refs if parent_id not in ids)
    return {
        "genes": genes,
        "mrnas": mrnas,
        "children": children,
        "feature_counts": feature_counts,
        "duplicate_ids": duplicate_ids,
        "missing_parent_refs": missing_parent_refs,
        "bad_coords": bad_coords,
    }


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            if end > merged[-1][1]:
                merged[-1][1] = end
    return [(start, end) for start, end in merged]


def load_te_intervals(path):
    intervals = defaultdict(list)
    if not path:
        return intervals
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return intervals
    with p.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            if start <= end:
                intervals[fields[0]].append((start, end))
    return {seqid: merge_intervals(vals) for seqid, vals in intervals.items()}


def overlap_bp(interval, intervals):
    if not intervals:
        return 0
    start, end = interval
    total = 0
    for te_start, te_end in intervals:
        if te_end < start:
            continue
        if te_start > end:
            break
        ov_start = max(start, te_start)
        ov_end = min(end, te_end)
        if ov_start <= ov_end:
            total += ov_end - ov_start + 1
    return total


def write_fasta(records, path):
    with Path(path).open("w", encoding="utf-8") as out:
        for name, desc, seq in records:
            out.write(">{} {}\n".format(name, desc).rstrip() + "\n")
            for line in wrap(seq):
                out.write(line + "\n")


def copy_nonempty_lines(src, dst):
    with Path(src).open(encoding="utf-8", errors="replace") as handle, \
            Path(dst).open("w", encoding="utf-8") as out:
        for line in handle:
            out.write(line)


def summarize_values(values):
    vals = sorted(values)
    if not vals:
        return {"count": 0, "min": 0, "median": 0, "mean": 0, "max": 0}
    n = len(vals)
    median = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2.0
    return {
        "count": n,
        "min": vals[0],
        "median": median,
        "mean": sum(vals) / n,
        "max": vals[-1],
    }


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome = read_fasta(args.genome)
    parsed = parse_gff3(args.gff3)
    te_intervals = load_te_intervals(args.te_gff3)
    genes = parsed["genes"]
    mrnas = parsed["mrnas"]
    children = parsed["children"]

    gene_gff = out_dir / "{}.genes.gff3".format(args.prefix)
    cds_fa = out_dir / "{}.cds.fa".format(args.prefix)
    protein_fa = out_dir / "{}.protein.fa".format(args.prefix)
    mrna_fa = out_dir / "{}.mrna.fa".format(args.prefix)
    qc_tsv = out_dir / "{}.gene_model_qc.tsv".format(args.prefix)
    low_conf_tsv = out_dir / "{}.low_confidence_gene_models.tsv".format(args.prefix)
    summary_tsv = out_dir / "{}.postprocess_summary.tsv".format(args.prefix)
    chr_tsv = out_dir / "{}.genes_by_sequence.qc.tsv".format(args.prefix)

    copy_nonempty_lines(args.gff3, gene_gff)

    cds_records = []
    protein_records = []
    mrna_records = []
    qc_rows = []
    low_conf_rows = []
    chr_gene_counts = Counter()
    flag_counts = Counter()
    gene_lengths = []
    cds_lengths = []
    protein_lengths = []
    exon_counts = []

    for mrna_id in sorted(mrnas):
        mrna = mrnas[mrna_id]
        parent_gene = mrna["attrs"].get("Parent", "")
        gene = genes.get(parent_gene, mrna)
        seqid = mrna["seqid"]
        strand = mrna["strand"]
        exons = children[mrna_id].get("exon", [])
        cdss = children[mrna_id].get("CDS", [])
        exon_intervals = [(rec["start"], rec["end"]) for rec in exons]
        cds_intervals = [(rec["start"], rec["end"]) for rec in cdss]
        gene_len = gene["end"] - gene["start"] + 1
        mrna_seq = extract_parts(genome, seqid, exon_intervals, strand)
        cds_seq = extract_parts(genome, seqid, cds_intervals, strand)
        protein_seq = translate(cds_seq)
        gene_te_bp = overlap_bp((gene["start"], gene["end"]), te_intervals.get(seqid, []))
        cds_te_bp = sum(overlap_bp(interval, te_intervals.get(seqid, [])) for interval in cds_intervals)
        cds_len = len(cds_seq)
        protein_len = len(protein_seq.rstrip("*"))
        internal_stop = "*" in protein_seq[:-1]
        has_stop = protein_seq.endswith("*")
        has_start = protein_seq.startswith("M")
        partial_codon = cds_len % 3 != 0
        gene_te_frac = gene_te_bp / gene_len if gene_len else 0
        cds_te_frac = cds_te_bp / cds_len if cds_len else 0
        flags = []
        if cds_len < args.short_cds:
            flags.append("short_cds")
        if len(exons) == 1:
            flags.append("single_exon")
        if partial_codon:
            flags.append("cds_not_multiple_of_3")
        if internal_stop:
            flags.append("internal_stop")
        if not has_start:
            flags.append("no_start_M")
        if not has_stop:
            flags.append("no_terminal_stop")
        if gene_te_frac >= args.gene_te_fraction:
            flags.append("high_gene_te_overlap")
        if cds_te_frac >= args.cds_te_fraction:
            flags.append("high_cds_te_overlap")
        if not flags:
            flags.append("PASS")
        for flag in flags:
            flag_counts[flag] += 1
        desc = "gene={} seqid={} strand={} cds_len={} flags={}".format(
            parent_gene, seqid, strand, cds_len, ",".join(flags)
        )
        if mrna_seq:
            mrna_records.append((mrna_id, desc, mrna_seq))
        if cds_seq:
            cds_records.append((mrna_id, desc, cds_seq))
            protein_records.append((mrna_id, desc, protein_seq.rstrip("*")))
        qc_rows.append([
            parent_gene,
            mrna_id,
            seqid,
            gene["start"],
            gene["end"],
            strand,
            gene_len,
            len(exons),
            len(cdss),
            len(mrna_seq),
            cds_len,
            protein_len,
            int(has_start),
            int(has_stop),
            int(internal_stop),
            int(partial_codon),
            gene_te_bp,
            "{:.6f}".format(gene_te_frac),
            cds_te_bp,
            "{:.6f}".format(cds_te_frac),
            ",".join(flags),
        ])
        if flags != ["PASS"]:
            low_conf_rows.append(qc_rows[-1])
        chr_gene_counts[seqid] += 1
        gene_lengths.append(gene_len)
        cds_lengths.append(cds_len)
        protein_lengths.append(protein_len)
        exon_counts.append(len(exons))

    write_fasta(mrna_records, mrna_fa)
    write_fasta(cds_records, cds_fa)
    write_fasta(protein_records, protein_fa)

    qc_header = [
        "gene_id", "mrna_id", "seqid", "start", "end", "strand", "gene_len",
        "exon_count", "cds_part_count", "mrna_len", "cds_len", "protein_len",
        "has_start_M", "has_terminal_stop", "internal_stop", "cds_not_multiple_of_3",
        "gene_te_overlap_bp", "gene_te_overlap_fraction",
        "cds_te_overlap_bp", "cds_te_overlap_fraction", "flags",
    ]
    with qc_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(qc_header)
        writer.writerows(qc_rows)
    with low_conf_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(qc_header)
        writer.writerows(low_conf_rows)
    with chr_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["seqid", "gene_count"])
        for seqid in sorted(chr_gene_counts):
            writer.writerow([seqid, chr_gene_counts[seqid]])

    summaries = {
        "gene_len": summarize_values(gene_lengths),
        "cds_len": summarize_values(cds_lengths),
        "protein_len": summarize_values(protein_lengths),
        "exon_count": summarize_values(exon_counts),
    }
    with summary_tsv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["genome_sequences", len(genome)])
        writer.writerow(["gff3_feature_gene", parsed["feature_counts"].get("gene", 0)])
        writer.writerow(["gff3_feature_mRNA", parsed["feature_counts"].get("mRNA", 0)])
        writer.writerow(["gff3_feature_exon", parsed["feature_counts"].get("exon", 0)])
        writer.writerow(["gff3_feature_CDS", parsed["feature_counts"].get("CDS", 0)])
        writer.writerow(["mrna_models", len(mrnas)])
        writer.writerow(["cds_fasta_records", len(cds_records)])
        writer.writerow(["protein_fasta_records", len(protein_records)])
        writer.writerow(["duplicate_ids", parsed["duplicate_ids"]])
        writer.writerow(["missing_parent_refs", parsed["missing_parent_refs"]])
        writer.writerow(["bad_coords", parsed["bad_coords"]])
        writer.writerow(["low_confidence_models", len(low_conf_rows)])
        for flag, count in sorted(flag_counts.items()):
            writer.writerow(["flag_{}".format(flag), count])
        for name, values in summaries.items():
            for key, value in values.items():
                writer.writerow(["{}_{}".format(name, key), value])

    print("Helixer postprocess/QC completed")
    print("GFF3:", gene_gff)
    print("CDS:", cds_fa)
    print("Protein:", protein_fa)
    print("mRNA:", mrna_fa)
    print("QC:", qc_tsv)
    print("Summary:", summary_tsv)


def parse_args():
    parser = argparse.ArgumentParser(description="Postprocess Helixer GFF3 and run basic QC.")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--te-gff3", default="")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--short-cds", type=int, default=300)
    parser.add_argument("--gene-te-fraction", type=float, default=0.8)
    parser.add_argument("--cds-te-fraction", type=float, default=0.5)
    return parser.parse_args()


if __name__ == "__main__":
    main()
