#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    name = None
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, header, "".join(chunks)


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def read_lengths(fasta):
    return {name: len(seq) for name, header, seq in read_fasta(fasta)}


def read_records(fasta):
    return {name: seq for name, header, seq in read_fasta(fasta)}


def parse_paf(path, source_label, min_mapq, min_aln_len):
    rows = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])
            if mapq < min_mapq or aln_len < min_aln_len:
                continue
            rows.append(
                {
                    "source": source_label,
                    "query": query,
                    "query_len": query_len,
                    "query_start": query_start,
                    "query_end": query_end,
                    "strand": strand,
                    "target": target,
                    "target_len": target_len,
                    "target_start": target_start,
                    "target_end": target_end,
                    "matches": matches,
                    "aln_len": aln_len,
                    "mapq": mapq,
                }
            )
    return rows


def interval_union_bp(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def load_pairs(path):
    pairs = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pairs.append(row)
    return pairs


def infer_assignments(rows, chr_lengths, pairs, dominance_threshold):
    intervals = defaultdict(lambda: {"hap1": [], "hap2": []})
    raw_bp = defaultdict(lambda: {"hap1": 0, "hap2": 0})
    best_targets = defaultdict(lambda: {"hap1": defaultdict(int), "hap2": defaultdict(int)})

    for row in rows:
        chrom = row["target"]
        if chrom not in chr_lengths:
            continue
        hap = row["source"]
        intervals[chrom][hap].append((row["target_start"], row["target_end"]))
        raw_bp[chrom][hap] += row["aln_len"]
        best_targets[chrom][hap][row["query"]] += row["aln_len"]

    summary = {}
    for chrom in sorted(chr_lengths):
        hap1_bp = interval_union_bp(intervals[chrom]["hap1"])
        hap2_bp = interval_union_bp(intervals[chrom]["hap2"])
        total = hap1_bp + hap2_bp
        if total == 0:
            dominant = "unknown"
            dominance = 0.0
        elif hap1_bp >= hap2_bp:
            dominant = "hap1"
            dominance = hap1_bp / total
        else:
            dominant = "hap2"
            dominance = hap2_bp / total
        status = "confident" if dominant != "unknown" and dominance >= dominance_threshold else "ambiguous"
        summary[chrom] = {
            "chr": chrom,
            "length": chr_lengths[chrom],
            "hap1_covered_bp": hap1_bp,
            "hap2_covered_bp": hap2_bp,
            "hap1_raw_aln_bp": raw_bp[chrom]["hap1"],
            "hap2_raw_aln_bp": raw_bp[chrom]["hap2"],
            "dominant_hap": dominant,
            "dominance": dominance,
            "status": status,
            "top_hap1_query": top_target(best_targets[chrom]["hap1"]),
            "top_hap2_query": top_target(best_targets[chrom]["hap2"]),
        }

    corrected_pairs = []
    for pair in pairs:
        a = pair["copy1_chr"]
        b = pair["copy2_chr"]
        a_hap = summary[a]["dominant_hap"]
        b_hap = summary[b]["dominant_hap"]

        if a_hap == "hap1" and b_hap == "hap2":
            copy1, copy2 = a, b
            action = "kept"
        elif a_hap == "hap2" and b_hap == "hap1":
            copy1, copy2 = b, a
            action = "swapped"
        else:
            copy1, copy2 = a, b
            action = "ambiguous_keep_original"

        corrected_pairs.append(
            {
                "pair_id": pair["pair_id"],
                "copy1_chr": copy1,
                "copy2_chr": copy2,
                "original_copy1_chr": a,
                "original_copy2_chr": b,
                "copy1_dominant_hap": summary[copy1]["dominant_hap"],
                "copy2_dominant_hap": summary[copy2]["dominant_hap"],
                "action": action,
            }
        )

    return summary, corrected_pairs


def top_target(counts):
    if not counts:
        return ""
    name, bp = max(counts.items(), key=lambda item: item[1])
    return f"{name}:{bp}"


def write_summary(path, summary):
    fields = [
        "chr",
        "length",
        "hap1_covered_bp",
        "hap2_covered_bp",
        "hap1_raw_aln_bp",
        "hap2_raw_aln_bp",
        "dominant_hap",
        "dominance",
        "status",
        "top_hap1_query",
        "top_hap2_query",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for chrom in sorted(summary):
            row = dict(summary[chrom])
            row["dominance"] = f"{row['dominance']:.6f}"
            writer.writerow(row)


def write_corrected_pairs(path, corrected_pairs):
    fields = [
        "pair_id",
        "copy1_chr",
        "copy2_chr",
        "original_copy1_chr",
        "original_copy2_chr",
        "copy1_dominant_hap",
        "copy2_dominant_hap",
        "action",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in corrected_pairs:
            writer.writerow(row)


def write_copy_fastas(records, corrected_pairs, copy1_out, copy2_out):
    with open(copy1_out, "w") as out1, open(copy2_out, "w") as out2:
        for idx, pair in enumerate(corrected_pairs, start=1):
            out_name = f"chr{idx:02d}"
            for out, source in [(out1, pair["copy1_chr"]), (out2, pair["copy2_chr"])]:
                out.write(f">{out_name}\n")
                for line in wrap(records[source]):
                    out.write(line + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Assign final chromosomes to hifiasm hap1/hap2 provenance.")
    parser.add_argument("--final-fasta", required=True)
    parser.add_argument("--hap1-paf", required=True)
    parser.add_argument("--hap2-paf", required=True)
    parser.add_argument("--homolog-pairs", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--copy1-fasta", required=True)
    parser.add_argument("--copy2-fasta", required=True)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--dominance-threshold", type=float, default=0.60)
    return parser.parse_args()


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    chr_lengths = read_lengths(args.final_fasta)
    records = read_records(args.final_fasta)
    pairs = load_pairs(args.homolog_pairs)
    rows = []
    rows.extend(parse_paf(args.hap1_paf, "hap1", args.min_mapq, args.min_aln_len))
    rows.extend(parse_paf(args.hap2_paf, "hap2", args.min_mapq, args.min_aln_len))

    summary, corrected_pairs = infer_assignments(
        rows,
        chr_lengths,
        pairs,
        args.dominance_threshold,
    )
    write_summary(str(out_prefix) + ".chr_hap_provenance.tsv", summary)
    write_corrected_pairs(str(out_prefix) + ".corrected_homolog_pairs.tsv", corrected_pairs)
    write_copy_fastas(records, corrected_pairs, args.copy1_fasta, args.copy2_fasta)

    ambiguous = [row for row in summary.values() if row["status"] != "confident"]
    swapped = [row for row in corrected_pairs if row["action"] == "swapped"]
    print(f"Wrote {out_prefix}.chr_hap_provenance.tsv")
    print(f"Wrote {out_prefix}.corrected_homolog_pairs.tsv")
    print(f"Wrote {args.copy1_fasta}")
    print(f"Wrote {args.copy2_fasta}")
    print(f"Ambiguous chromosomes: {len(ambiguous)}")
    print(f"Swapped homolog pairs: {len(swapped)}")


if __name__ == "__main__":
    main()
