#!/usr/bin/env python3
"""Structure and splice-support QC for a genome-anchored GFF3 annotation.

Introns are derived from adjacent exons of each transcript, as in the
reference methods.  The program intentionally does not infer introns from
Trinity contigs, because those contigs do not carry reliable genome
coordinates.
"""

import argparse
import bisect
import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


TRANSCRIPT_FEATURES = {"mrna", "transcript"}
EXON_FEATURES = {"exon"}
CDS_FEATURES = {"cds"}
REF_CONSUMING = {"M", "D", "N", "=", "X"}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_attributes(text):
    result = {}
    for raw in text.strip().rstrip(";").split(";"):
        item = raw.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
            value = value.strip().strip('"')
        result[key] = value
    return result


def parse_gff3(path):
    genes = {}
    transcripts = {}
    exons_by_tx = defaultdict(list)
    cds_by_tx = defaultdict(list)
    feature_counts = Counter()
    feature_bp = Counter()
    feature_lengths = defaultdict(list)
    all_ids = set()
    duplicate_ids = 0
    bad_rows = 0
    bad_coords = 0
    parent_refs = set()

    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                bad_rows += 1
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs_text = fields
            feature_key = feature.lower()
            feature_counts[feature] += 1
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                bad_coords += 1
                continue
            if start_i < 1 or end_i < start_i:
                bad_coords += 1
                continue
            feature_bp[feature] += end_i - start_i + 1
            feature_lengths[feature_key].append(end_i - start_i + 1)
            attrs = parse_attributes(attrs_text)
            record = {
                "seqid": seqid,
                "start": start_i,
                "end": end_i,
                "strand": strand if strand in {"+", "-"} else ".",
                "feature": feature,
                "attrs": attrs,
                "line": line_number,
            }
            record_id = attrs.get("ID")
            if record_id:
                if record_id in all_ids:
                    duplicate_ids += 1
                all_ids.add(record_id)
            if feature_key == "gene":
                if record_id:
                    genes[record_id] = record
            elif feature_key in TRANSCRIPT_FEATURES:
                tx_id = record_id or attrs.get("transcript_id")
                if tx_id:
                    record["id"] = tx_id
                    record["gene_id"] = attrs.get("Parent", attrs.get("gene_id", "")).split(",")[0]
                    transcripts[tx_id] = record
            elif feature_key in EXON_FEATURES or feature_key in CDS_FEATURES:
                parents = attrs.get("Parent", attrs.get("transcript_id", ""))
                for parent in filter(None, (p.strip() for p in parents.split(","))):
                    parent_refs.add(parent)
                    if feature_key in EXON_FEATURES:
                        exons_by_tx[parent].append(record)
                    else:
                        cds_by_tx[parent].append(record)

    for tx_id, tx in transcripts.items():
        gene_id = tx.get("gene_id", "")
        if gene_id and gene_id not in genes:
            # Some callers use a gene_id attribute without a separate gene row.
            genes.setdefault(gene_id, {
                "seqid": tx["seqid"], "start": tx["start"], "end": tx["end"],
                "strand": tx["strand"], "feature": "gene", "attrs": {}, "line": tx["line"],
            })
    return {
        "genes": genes,
        "transcripts": transcripts,
        "exons": exons_by_tx,
        "cds": cds_by_tx,
        "feature_counts": feature_counts,
        "feature_bp": feature_bp,
        "feature_lengths": feature_lengths,
        "duplicate_ids": duplicate_ids,
        "bad_rows": bad_rows,
        "bad_coords": bad_coords,
        "parent_refs": parent_refs,
    }


def derive_structures(parsed):
    transcript_rows = []
    exon_rows = []
    intron_rows = []
    missing_parent = 0
    exon_overlap = 0
    invalid_introns = 0
    gene_tx_counts = Counter()
    intron_key_to_ids = defaultdict(list)

    for tx_id, tx in sorted(parsed["transcripts"].items()):
        gene_id = tx.get("gene_id", "")
        exons = list(parsed["exons"].get(tx_id, []))
        cds = list(parsed["cds"].get(tx_id, []))
        if not exons:
            missing_parent += 1
        strand = tx["strand"]
        oriented = sorted(exons, key=lambda r: (r["start"], r["end"]), reverse=(strand == "-"))
        gene_tx_counts[gene_id] += 1
        tx_start = min([tx["start"]] + [r["start"] for r in exons])
        tx_end = max([tx["end"]] + [r["end"] for r in exons])
        exon_bp = sum(r["end"] - r["start"] + 1 for r in exons)
        cds_bp = sum(r["end"] - r["start"] + 1 for r in cds)
        transcript_rows.append({
            "transcript_id": tx_id,
            "gene_id": gene_id,
            "seqid": tx["seqid"],
            "strand": strand,
            "start": tx_start,
            "end": tx_end,
            "transcript_span": tx_end - tx_start + 1,
            "exon_count": len(exons),
            "intron_count": max(0, len(exons) - 1),
            "exon_bp": exon_bp,
            "cds_bp": cds_bp,
            "cds_mod3": cds_bp % 3,
            "has_cds": int(bool(cds)),
        })

        for rank, exon in enumerate(oriented, 1):
            exon_id = exon["attrs"].get("ID", f"{tx_id}.exon{rank}")
            exon_rows.append({
                "exon_id": exon_id,
                "transcript_id": tx_id,
                "gene_id": gene_id,
                "seqid": exon["seqid"],
                "strand": strand,
                "rank": rank,
                "start": exon["start"],
                "end": exon["end"],
                "length": exon["end"] - exon["start"] + 1,
            })
            if rank > 1:
                previous = oriented[rank - 2]
                if strand == "-":
                    intron_start = oriented[rank - 1]["end"] + 1
                    intron_end = previous["start"] - 1
                else:
                    intron_start = previous["end"] + 1
                    intron_end = oriented[rank - 1]["start"] - 1
                if intron_end < intron_start:
                    invalid_introns += 1
                    continue
                if (strand != "-" and previous["end"] >= oriented[rank - 1]["start"]) or (
                    strand == "-" and oriented[rank - 1]["end"] >= previous["start"]
                ):
                    exon_overlap += 1
                intron_id = f"{tx_id}.intron{rank - 1}"
                row = {
                    "intron_id": intron_id,
                    "transcript_id": tx_id,
                    "gene_id": gene_id,
                    "seqid": tx["seqid"],
                    "strand": strand,
                    "rank": rank - 1,
                    "start": min(intron_start, intron_end),
                    "end": max(intron_start, intron_end),
                    "length": abs(intron_end - intron_start) + 1,
                    "previous_exon_id": previous["attrs"].get("ID", ""),
                    "next_exon_id": oriented[rank - 1]["attrs"].get("ID", ""),
                }
                intron_rows.append(row)
                intron_key_to_ids[(row["seqid"], row["start"], row["end"])].append(intron_id)

    return {
        "transcripts": transcript_rows,
        "exons": exon_rows,
        "introns": intron_rows,
        "intron_key_to_ids": intron_key_to_ids,
        "missing_parent": missing_parent,
        "exon_overlap": exon_overlap,
        "invalid_introns": invalid_introns,
        "gene_tx_counts": gene_tx_counts,
    }


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_bam_specs(specs):
    result = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"--bam must use LABEL=PATH: {spec}")
        label, path = spec.split("=", 1)
        if not label or not path:
            raise ValueError(f"--bam must use LABEL=PATH: {spec}")
        result.append((label, Path(path)))
    return result


def iter_sam_records(bam, threads):
    # Exclude unmapped, secondary and supplementary alignments so one read
    # cannot inflate splice support through alternate alignment records.
    command = ["samtools", "view", "-@", str(max(1, threads)), "-F", "2308", str(bam)]
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except OSError as exc:
        raise RuntimeError(f"could not start samtools for {bam}: {exc}") from exc
    assert process.stdout is not None
    for line in process.stdout:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 6 and fields[2] != "*" and fields[5] != "*":
            yield fields[0], fields[2], int(fields[3]), fields[5]
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"samtools view failed for {bam}: {stderr.strip()}")


def splice_support(bam_specs, intron_keys, threads):
    support = {label: Counter() for label, _ in bam_specs}
    for label, bam in bam_specs:
        for qname, seqid, pos, cigar in iter_sam_records(bam, threads):
            ref_pos = pos
            for length_text, op in CIGAR_RE.findall(cigar):
                length = int(length_text)
                if op == "N":
                    key = (seqid, ref_pos, ref_pos + length - 1)
                    if key in intron_keys:
                        support[label][key] += 1
                    ref_pos += length
                elif op in REF_CONSUMING:
                    ref_pos += length
    rows = []
    for key in sorted(intron_keys):
        seqid, start, end = key
        row = {"seqid": seqid, "start": start, "end": end, "length": end - start + 1,
               "transcript_count": len(intron_keys[key])}
        for label, _ in bam_specs:
            row[f"{label}_supporting_alignments"] = support[label].get(key, 0)
        rows.append(row)
    return rows


def load_te_intervals(path):
    intervals = defaultdict(list)
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if end >= start:
                intervals[fields[0]].append((start, end, fields[2]))
    for seqid in intervals:
        intervals[seqid].sort()
    return intervals


def te_overlap(introns, te_path):
    intervals = load_te_intervals(te_path)
    by_seq = {}
    for seqid, values in intervals.items():
        prefix_max_end = []
        current_max = 0
        for _, end, _ in values:
            current_max = max(current_max, end)
            prefix_max_end.append(current_max)
        by_seq[seqid] = (values, [item[0] for item in values], prefix_max_end)
    rows = []
    for intron in introns:
        values, starts, prefix_max_end = by_seq.get(intron["seqid"], ([], [], []))
        # The prefix maximum lets us skip intervals that end before the
        # intron, while retaining intervals that start upstream but overlap it.
        left = bisect.bisect_left(prefix_max_end, intron["start"])
        right = bisect.bisect_right(starts, intron["end"])
        overlaps = []
        te_count = 0
        te_types = set()
        for idx in range(left, right):
            ov_start = max(intron["start"], values[idx][0])
            ov_end = min(intron["end"], values[idx][1])
            ov = ov_end - ov_start + 1
            if ov > 0:
                overlaps.append((ov_start, ov_end))
                te_count += 1
                te_types.add(values[idx][2])
        # EDTA features can overlap; report the union rather than double
        # counting bases covered by multiple TE records.
        overlap_bp = 0
        for ov_start, ov_end in sorted(overlaps):
            if overlap_bp == 0:
                merged_start, merged_end = ov_start, ov_end
                overlap_bp = ov_end - ov_start + 1
            elif ov_start <= merged_end + 1:
                if ov_end > merged_end:
                    overlap_bp += ov_end - merged_end
                    merged_end = ov_end
            else:
                merged_start, merged_end = ov_start, ov_end
                overlap_bp += ov_end - ov_start + 1
        rows.append({
            "intron_id": intron["intron_id"], "transcript_id": intron["transcript_id"],
            "gene_id": intron["gene_id"], "seqid": intron["seqid"],
            "start": intron["start"], "end": intron["end"], "length": intron["length"],
            "te_overlap_bp": overlap_bp,
            "te_overlap_fraction": f"{overlap_bp / intron['length']:.8f}",
            "te_feature_count": te_count, "te_feature_types": ",".join(sorted(te_types)),
        })
    return rows


def build_summary(parsed, structures, bam_specs, te_rows):
    counts = parsed["feature_counts"]
    bp = parsed["feature_bp"]
    tx_rows = structures["transcripts"]
    introns = structures["introns"]
    exon_count = sum(1 for row in structures["exons"])
    cds_count = sum(count for feature, count in counts.items() if feature.lower() in CDS_FEATURES)
    rna_number = len(parsed["transcripts"])
    long_1k = sum(row["length"] >= 1000 for row in introns)
    long_10k = sum(row["length"] >= 10000 for row in introns)
    te_overlap_count = sum(int(row["te_overlap_bp"]) > 0 for row in te_rows) if te_rows is not None else 0
    exon_lengths = parsed["feature_lengths"].get("exon", [])
    cds_lengths = parsed["feature_lengths"].get("cds", [])
    rows = [
        {"metric": "Gene number", "value": len(parsed["genes"]), "notes": "unique gene IDs"},
        {"metric": "RNA number", "value": rna_number, "notes": "unique mRNA/transcript IDs"},
        {"metric": "CDS number", "value": cds_count, "notes": "CDS feature rows"},
        {"metric": "Exon number", "value": exon_count, "notes": "exon feature rows"},
        {"metric": "Intron number", "value": len(introns), "notes": "derived from adjacent exons"},
        {"metric": "Genes with transcript", "value": sum(bool(gene_id) and bool(number) for gene_id, number in structures["gene_tx_counts"].items()), "notes": "genes represented by transcript records"},
        {"metric": "Mean exons per transcript", "value": f"{sum(r['exon_count'] for r in tx_rows) / len(tx_rows):.4f}" if tx_rows else "0", "notes": "derived"},
        {"metric": "Mean CDS length (bp)", "value": f"{sum(r['cds_bp'] for r in tx_rows if r['cds_bp']) / max(1, sum(r['cds_bp'] > 0 for r in tx_rows)):.2f}", "notes": "per transcript with CDS"},
        {"metric": "Mean exon length (bp)", "value": f"{sum(exon_lengths) / len(exon_lengths):.2f}" if exon_lengths else "0", "notes": "per exon feature"},
        {"metric": "Mean intron length (bp)", "value": f"{sum(row['length'] for row in introns) / len(introns):.2f}" if introns else "0", "notes": "derived per intron"},
        {"metric": "Long introns >=1 kb", "value": long_1k, "notes": "candidate long introns"},
        {"metric": "Long introns >=10 kb", "value": long_10k, "notes": "candidate very long introns"},
        {"metric": "TE-overlapping introns", "value": te_overlap_count, "notes": "only when --te-gff3 supplied"},
        {"metric": "GFF3 duplicate IDs", "value": parsed["duplicate_ids"], "notes": "duplicate ID attributes"},
        {"metric": "GFF3 malformed rows", "value": parsed["bad_rows"], "notes": "rows not containing 9 columns"},
        {"metric": "GFF3 invalid coordinates", "value": parsed["bad_coords"], "notes": "non-integer or end < start"},
        {"metric": "Missing exon parents", "value": structures["missing_parent"], "notes": "transcripts without exon records"},
        {"metric": "Overlapping adjacent exons", "value": structures["exon_overlap"], "notes": "coordinate QC flag"},
        {"metric": "Invalid derived introns", "value": structures["invalid_introns"], "notes": "negative or zero-length introns"},
    ]
    return rows


def parse_args():
    parser = argparse.ArgumentParser(description="Derive introns and summarize GFF3 gene structures.")
    parser.add_argument("--gff3", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--bam", action="append", default=[], metavar="LABEL=PATH", help="repeat for splice support")
    parser.add_argument("--te-gff3", type=Path, help="optional EDTA TE GFF3 for intron overlap")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--long-intron-thresholds", default="1000,10000")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.gff3.is_file() or args.gff3.stat().st_size == 0:
        raise SystemExit(f"ERROR: GFF3 missing or empty: {args.gff3}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    try:
        bam_specs = parse_bam_specs(args.bam)
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}")
    for label, path in bam_specs:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: BAM missing or empty for {label}: {path}")
    if args.te_gff3 and (not args.te_gff3.is_file() or args.te_gff3.stat().st_size == 0):
        raise SystemExit(f"ERROR: TE GFF3 missing or empty: {args.te_gff3}")

    parsed = parse_gff3(args.gff3)
    structures = derive_structures(parsed)
    intron_keys = structures["intron_key_to_ids"]
    te_rows = te_overlap(structures["introns"], args.te_gff3) if args.te_gff3 else None

    tx_cols = ["transcript_id", "gene_id", "seqid", "strand", "start", "end", "transcript_span", "exon_count", "intron_count", "exon_bp", "cds_bp", "cds_mod3", "has_cds"]
    exon_cols = ["exon_id", "transcript_id", "gene_id", "seqid", "strand", "rank", "start", "end", "length"]
    intron_cols = ["intron_id", "transcript_id", "gene_id", "seqid", "strand", "rank", "start", "end", "length", "previous_exon_id", "next_exon_id"]
    write_tsv(args.out_dir / "transcript_structure.tsv", tx_cols, structures["transcripts"])
    write_tsv(args.out_dir / "exon_details.tsv", exon_cols, structures["exons"])
    write_tsv(args.out_dir / "intron_details.tsv", intron_cols, structures["introns"])

    thresholds = []
    for token in args.long_intron_thresholds.split(","):
        if token.strip().isdigit():
            thresholds.append(int(token.strip()))
    thresholds = sorted(set(thresholds)) or [1000, 10000]
    long_rows = []
    for threshold in thresholds:
        for row in structures["introns"]:
            if row["length"] >= threshold:
                long_rows.append({"threshold_bp": threshold, **row})
    write_tsv(args.out_dir / "long_intron_candidates.tsv", ["threshold_bp"] + intron_cols, long_rows)

    length_bins = [(0, 99, "0-99"), (100, 499, "100-499"), (500, 999, "500-999"), (1000, 4999, "1000-4999"), (5000, 9999, "5000-9999"), (10000, 49999, "10000-49999"), (50000, None, "50000+")]
    length_rows = []
    total_introns = len(structures["introns"])
    for low, high, label in length_bins:
        count = sum(low <= row["length"] <= high if high is not None else row["length"] >= low for row in structures["introns"])
        length_rows.append({"length_bin": label, "lower_bp": low, "upper_bp": high if high is not None else "", "intron_count": count, "fraction": f"{count / total_introns:.8f}" if total_introns else "0"})
    write_tsv(args.out_dir / "intron_length_distribution.tsv", ["length_bin", "lower_bp", "upper_bp", "intron_count", "fraction"], length_rows)

    def write_feature_distribution(name, lengths):
        total = len(lengths)
        rows = []
        for low, high, label in length_bins:
            count = sum(low <= length <= high if high is not None else length >= low for length in lengths)
            rows.append({"length_bin": label, "lower_bp": low, "upper_bp": high if high is not None else "", "feature_count": count, "fraction": f"{count / total:.8f}" if total else "0"})
        write_tsv(args.out_dir / name, ["length_bin", "lower_bp", "upper_bp", "feature_count", "fraction"], rows)

    write_feature_distribution("exon_length_distribution.tsv", parsed["feature_lengths"].get("exon", []))
    write_feature_distribution("cds_length_distribution.tsv", parsed["feature_lengths"].get("cds", []))
    exon_number_counts = Counter(row["exon_count"] for row in structures["transcripts"])
    write_tsv(args.out_dir / "transcript_exon_count_distribution.tsv", ["exon_count", "transcript_count", "fraction"], [{"exon_count": count, "transcript_count": number, "fraction": f"{number / len(structures['transcripts']):.8f}" if structures["transcripts"] else "0"} for count, number in sorted(exon_number_counts.items())])

    support_rows = None
    if bam_specs:
        support_rows = splice_support(bam_specs, intron_keys, args.threads)
        write_tsv(args.out_dir / "splice_junction_support.tsv", ["seqid", "start", "end", "length", "transcript_count"] + [f"{label}_supporting_alignments" for label, _ in bam_specs], support_rows)
    else:
        write_tsv(args.out_dir / "splice_junction_support.tsv", ["seqid", "start", "end", "length", "transcript_count"], [{"seqid": k[0], "start": k[1], "end": k[2], "length": k[2] - k[1] + 1, "transcript_count": len(v)} for k, v in sorted(intron_keys.items())])
    if te_rows is not None:
        write_tsv(args.out_dir / "intron_te_overlap.tsv", ["intron_id", "transcript_id", "gene_id", "seqid", "start", "end", "length", "te_overlap_bp", "te_overlap_fraction", "te_feature_count", "te_feature_types"], te_rows)

    summary = build_summary(parsed, structures, bam_specs, te_rows)
    if support_rows is not None:
        for label, _ in bam_specs:
            column = f"{label}_supporting_alignments"
            supported = sum(int(row[column]) > 0 for row in support_rows)
            summary.append({"metric": f"Introns with splice support ({label})", "value": supported, "notes": "exact annotated junction match"})
    write_tsv(args.out_dir / "gene_structure_summary.tsv", ["metric", "value", "notes"], summary)
    summary_values = {row["metric"]: row["value"] for row in summary}
    write_tsv(
        args.out_dir / "gene_structure_metrics.tsv",
        ["species", "Gene number", "RNA number", "CDS number", "Exon number", "Intron number", "Mean exons per transcript", "Mean CDS length (bp)", "Mean exon length (bp)", "Mean intron length (bp)", "Long introns >=1 kb", "Long introns >=10 kb", "TE-overlapping introns"],
        [{"species": args.prefix, **{key: summary_values.get(key, "") for key in ["Gene number", "RNA number", "CDS number", "Exon number", "Intron number", "Mean exons per transcript", "Mean CDS length (bp)", "Mean exon length (bp)", "Mean intron length (bp)", "Long introns >=1 kb", "Long introns >=10 kb", "TE-overlapping introns"]}}],
    )

    with (args.out_dir / "structure_qc_manifest.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"GFF3: {args.gff3}\n")
        handle.write(f"Prefix: {args.prefix}\n")
        handle.write(f"TE GFF3: {args.te_gff3 if args.te_gff3 else 'none'}\n")
        handle.write("BAM splice evidence:\n")
        for label, path in bam_specs:
            handle.write(f"  {label}\t{path}\n")
        handle.write(f"Long intron thresholds: {','.join(str(x) for x in thresholds)}\n")
        handle.write("Intron definition: adjacent exons within each transcript; intron_start=exon_i.end+1; intron_end=exon_(i+1).start-1.\n")
        handle.write("Splice definition: exact seqid/start/end match to an N CIGAR operation in samtools view output.\n")
    print(f"Structure QC completed: {args.out_dir}")


if __name__ == "__main__":
    main()
