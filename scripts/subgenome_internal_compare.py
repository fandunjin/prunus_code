#!/usr/bin/env python3
"""Summarize final x=8 haplotype/subgenome comparisons.

The script has two modes:
  paf-summary: convert a minimap2 PAF file into a chromosome-pair summary.
  collect: combine synteny summaries, BUSCO summaries, and RBH protein matches.
"""

import argparse
import csv
import re
import shutil
from pathlib import Path


HAPS = ["hapA1", "hapA2", "hapB1", "hapB2"]
GROUPS = {"A": ["hapA1", "hapA2"], "B": ["hapB1", "hapB2"]}
PAIRWISE_RBH = [
    ("hapA1", "hapA2"),
    ("hapB1", "hapB2"),
    ("hapA1", "hapB1"),
    ("hapA1", "hapB2"),
    ("hapA2", "hapB1"),
    ("hapA2", "hapB2"),
]
PACKAGE_FILES = {
    "genome": ("genome", ".fa"),
    "gff": ("gff", ".gff3"),
    "cds": ("cds", ".cds.fa"),
    "pep": ("pep", ".pep.fa"),
}


def natural_key(text):
    parts = []
    buf = ""
    digit = False
    for ch in text:
        if buf and ch.isdigit() != digit:
            parts.append(int(buf) if digit else buf)
            buf = ch
        else:
            buf += ch
        digit = ch.isdigit()
    if buf:
        parts.append(int(buf) if digit else buf)
    return parts


def read_tsv(path):
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def package_input(package_dir, hap, kind):
    """Resolve one haplotype input from the renamed or legacy flat package."""
    subdir, suffix = PACKAGE_FILES[kind]
    label = "pzi_%s" % hap
    candidates = [
        package_dir / subdir / (label + suffix),
        package_dir / (label + suffix),
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    raise FileNotFoundError(
        "No usable %s input for %s; checked: %s"
        % (kind, hap, ", ".join(str(path) for path in candidates))
    )


def read_fasta_ids_and_lengths(path):
    ids = []
    lengths = {}
    current = None
    total = 0
    if not path.is_file():
        return ids, lengths
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = total
                current = line[1:].strip().split()[0]
                ids.append(current)
                total = 0
            else:
                total += len(line.strip())
    if current is not None:
        lengths[current] = total
    return ids, lengths


def parse_gff_attributes(text):
    attrs = {}
    for part in text.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


def summarize_gff(path):
    counts = {"gene": 0, "mRNA": 0, "CDS": 0, "exon": 0}
    seq_counts = {}
    gene_ids = set()
    if not path.is_file():
        return counts, seq_counts, gene_ids
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, feature = fields[0], fields[2]
            seq_counts[seqid] = seq_counts.get(seqid, 0) + 1
            if feature in counts:
                counts[feature] += 1
            if feature == "gene":
                attrs = parse_gff_attributes(fields[8])
                gene_ids.add(attrs.get("ID", "gene_%s" % counts["gene"]))
    return counts, seq_counts, gene_ids


def parse_paf(path):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                rows.append({
                    "query": fields[0],
                    "query_len": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "strand": fields[4],
                    "target": fields[5],
                    "target_len": int(fields[6]),
                    "target_start": int(fields[7]),
                    "target_end": int(fields[8]),
                    "matches": int(fields[9]),
                    "block_len": int(fields[10]),
                    "mapq": int(fields[11]),
                })
            except ValueError:
                continue
    return rows


def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            if end > cur_end:
                cur_end = end
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def write_paf_summary(paf, out_path):
    rows = parse_paf(paf)
    grouped = {}
    target_lengths = {}
    query_lengths = {}
    for row in rows:
        key = (row["target"], row["query"])
        target_lengths[row["target"]] = row["target_len"]
        query_lengths[row["query"]] = row["query_len"]
        data = grouped.setdefault(key, {
            "alignment_count": 0,
            "aligned_block_bp": 0,
            "matching_bp": 0,
            "target_intervals": [],
            "query_intervals": [],
        })
        data["alignment_count"] += 1
        data["aligned_block_bp"] += row["block_len"]
        data["matching_bp"] += row["matches"]
        data["target_intervals"].append((row["target_start"], row["target_end"]))
        data["query_intervals"].append((row["query_start"], row["query_end"]))

    out_rows = []
    for target, query in sorted(grouped, key=lambda x: (natural_key(x[0]), natural_key(x[1]))):
        data = grouped[(target, query)]
        target_covered = merge_intervals(data["target_intervals"])
        query_covered = merge_intervals(data["query_intervals"])
        block_bp = data["aligned_block_bp"]
        identity = data["matching_bp"] / block_bp if block_bp else 0
        out_rows.append({
            "target_chr": target,
            "query_chr": query,
            "target_len": target_lengths.get(target, ""),
            "query_len": query_lengths.get(query, ""),
            "alignment_count": data["alignment_count"],
            "aligned_block_bp": block_bp,
            "matching_bp": data["matching_bp"],
            "mean_identity": "%.6f" % identity,
            "target_covered_bp": target_covered,
            "query_covered_bp": query_covered,
            "target_covered_pct": "%.4f" % (100 * target_covered / target_lengths[target] if target_lengths.get(target) else 0),
            "query_covered_pct": "%.4f" % (100 * query_covered / query_lengths[query] if query_lengths.get(query) else 0),
        })
    write_tsv(out_path, out_rows, [
        "target_chr", "query_chr", "target_len", "query_len", "alignment_count",
        "aligned_block_bp", "matching_bp", "mean_identity",
        "target_covered_bp", "query_covered_bp", "target_covered_pct", "query_covered_pct",
    ])


def find_busco_summaries(busco_dir):
    summaries = {}
    if not busco_dir.is_dir():
        return summaries
    for path in busco_dir.rglob("short_summary*.txt"):
        hap = None
        text_path = str(path)
        for candidate in HAPS:
            if candidate in text_path:
                hap = candidate
                break
        if hap:
            summaries[hap] = parse_busco(path)
    return summaries


def parse_busco(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    row = {
        "busco_path": str(path),
        "complete_pct": "",
        "single_pct": "",
        "duplicated_pct": "",
        "fragmented_pct": "",
        "missing_pct": "",
        "n_buscos": "",
        "status": "parsed_without_compact_line",
    }
    m = re.search(r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:(\d+)", text)
    if m:
        row.update({
            "complete_pct": m.group(1),
            "single_pct": m.group(2),
            "duplicated_pct": m.group(3),
            "fragmented_pct": m.group(4),
            "missing_pct": m.group(5),
            "n_buscos": m.group(6),
            "status": "pass",
        })
    return row


def read_best_hits(path):
    best = {}
    if not path.is_file():
        return best
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            qid, sid = fields[0], fields[1]
            try:
                pident = float(fields[2])
                length = int(float(fields[3]))
                evalue = fields[4]
                bitscore = float(fields[5])
            except ValueError:
                continue
            old = best.get(qid)
            if old is None or bitscore > old["bitscore"]:
                best[qid] = {
                    "subject": sid,
                    "pident": pident,
                    "length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                }
    return best


def rbh_for_pair(blast_dir, query_hap, target_hap, protein_ids):
    forward = read_best_hits(blast_dir / ("%s_vs_%s.blastp.tsv" % (query_hap, target_hap)))
    reverse = read_best_hits(blast_dir / ("%s_vs_%s.blastp.tsv" % (target_hap, query_hap)))
    rbh = []
    for qid, hit in forward.items():
        sid = hit["subject"]
        rev = reverse.get(sid)
        if rev and rev["subject"] == qid:
            rbh.append({
                "query_gene": qid,
                "target_gene": sid,
                "pident": "%.3f" % hit["pident"],
                "length": hit["length"],
                "evalue": hit["evalue"],
                "bitscore": "%.1f" % hit["bitscore"],
            })
    query_ids = set(protein_ids.get(query_hap, []))
    target_ids = set(protein_ids.get(target_hap, []))
    retained_query = set(row["query_gene"] for row in rbh)
    retained_target = set(row["target_gene"] for row in rbh)
    return rbh, sorted(query_ids - retained_query, key=natural_key), sorted(target_ids - retained_target, key=natural_key)


def hap_metric(package_dir, hap):
    fa = package_input(package_dir, hap, "genome")
    gff = package_input(package_dir, hap, "gff")
    cds = package_input(package_dir, hap, "cds")
    pep = package_input(package_dir, hap, "pep")
    genome_ids, genome_lengths = read_fasta_ids_and_lengths(fa)
    cds_ids, cds_lengths = read_fasta_ids_and_lengths(cds)
    pep_ids, pep_lengths = read_fasta_ids_and_lengths(pep)
    gff_counts, seq_counts, gene_ids = summarize_gff(gff)
    return {
        "hap_set": hap,
        "genome_sequences": len(genome_ids),
        "genome_bp": sum(genome_lengths.values()),
        "gff_gene": gff_counts["gene"],
        "gff_mRNA": gff_counts["mRNA"],
        "gff_CDS": gff_counts["CDS"],
        "gff_exon": gff_counts["exon"],
        "cds_records": len(cds_ids),
        "cds_bp": sum(cds_lengths.values()),
        "pep_records": len(pep_ids),
        "pep_aa": sum(pep_lengths.values()),
    }, pep_ids


def read_synteny_summaries(work_dir):
    rows = []
    for path in sorted((work_dir / "synteny").glob("*.summary.tsv")):
        summary_rows = read_tsv(path)
        total_block = sum(int(float(row.get("aligned_block_bp", 0) or 0)) for row in summary_rows)
        total_target_cov = sum(int(float(row.get("target_covered_bp", 0) or 0)) for row in summary_rows)
        total_query_cov = sum(int(float(row.get("query_covered_bp", 0) or 0)) for row in summary_rows)
        pairs = len(summary_rows)
        prefix = path.name.replace(".summary.tsv", "")
        rows.append({
            "comparison": prefix,
            "chromosome_pairs_with_alignment": pairs,
            "total_aligned_block_bp": total_block,
            "total_target_covered_bp": total_target_cov,
            "total_query_covered_bp": total_query_cov,
            "summary_path": str(path),
        })
    return rows


def score_hap_for_representative(hap, metrics, busco):
    b = busco.get(hap, {})
    if b.get("complete_pct") != "":
        c = float(b.get("complete_pct", 0) or 0)
        m = float(b.get("missing_pct", 0) or 0)
        d = float(b.get("duplicated_pct", 0) or 0)
        return c - m - 0.20 * d, "BUSCO complete - missing - 0.2*duplicated"
    metric = metrics[hap]
    return float(metric["pep_records"]), "protein gene count fallback; BUSCO missing"


def concatenate_files(out_path, input_paths, gff=False):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as out:
        if gff:
            out.write("##gff-version 3\n")
        for path in input_paths:
            with path.open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if gff and line.startswith("#"):
                        continue
                    out.write(line)
                    if not line.endswith("\n"):
                        out.write("\n")


def write_representative_package(package_dir, out_dir, chosen_a, chosen_b):
    rep_dir = out_dir / "downstream_representative_AplusB"
    rep_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "fa": [package_input(package_dir, hap, "genome") for hap in (chosen_a, chosen_b)],
        "gff3": [package_input(package_dir, hap, "gff") for hap in (chosen_a, chosen_b)],
        "cds.fa": [package_input(package_dir, hap, "cds") for hap in (chosen_a, chosen_b)],
        "pep.fa": [package_input(package_dir, hap, "pep") for hap in (chosen_a, chosen_b)],
    }
    outputs = {
        "fa": rep_dir / "pzi_representative_AplusB.fa",
        "gff3": rep_dir / "pzi_representative_AplusB.gff3",
        "cds.fa": rep_dir / "pzi_representative_AplusB.cds.fa",
        "pep.fa": rep_dir / "pzi_representative_AplusB.pep.fa",
    }
    for kind, paths in inputs.items():
        concatenate_files(outputs[kind], paths, gff=(kind == "gff3"))
    rows = []
    for kind, out_path in outputs.items():
        rows.append({
            "file_type": kind,
            "path": str(out_path),
            "source_haplotypes": "%s,%s" % (chosen_a, chosen_b),
            "source_files": ",".join(str(p) for p in inputs[kind]),
        })
    write_tsv(rep_dir / "MANIFEST.tsv", rows, ["file_type", "path", "source_haplotypes", "source_files"])
    (rep_dir / "README.md").write_text(
        "Representative A+B haplotype package for species-level comparative genomics.\n\n"
        "This package concatenates one A haplotype and one B haplotype, preserving the full inferred A+B gene complement.\n"
        "Do not use hapA1 alone as the species-level representative unless the analysis is explicitly A-subgenome-only.\n",
        encoding="utf-8",
    )


def collect(args):
    package_dir = args.package_dir
    out_dir = args.out_dir
    work_dir = args.work_dir
    blast_dir = work_dir / "protein_rbh"
    busco_dir = work_dir / "busco"
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = {}
    protein_ids = {}
    for hap in HAPS:
        metric, pep_ids = hap_metric(package_dir, hap)
        metrics[hap] = metric
        protein_ids[hap] = pep_ids
    write_tsv(out_dir / "hap_gene_count_summary.tsv", [metrics[h] for h in HAPS], [
        "hap_set", "genome_sequences", "genome_bp", "gff_gene", "gff_mRNA",
        "gff_CDS", "gff_exon", "cds_records", "cds_bp", "pep_records", "pep_aa",
    ])

    group_rows = []
    for group, haps in GROUPS.items():
        group_rows.append({
            "group": group,
            "hap_sets": ",".join(haps),
            "genome_bp": sum(metrics[h]["genome_bp"] for h in haps),
            "gff_gene": sum(metrics[h]["gff_gene"] for h in haps),
            "pep_records": sum(metrics[h]["pep_records"] for h in haps),
            "cds_bp": sum(metrics[h]["cds_bp"] for h in haps),
        })
    a_gene = group_rows[0]["gff_gene"]
    b_gene = group_rows[1]["gff_gene"]
    group_rows.append({
        "group": "A_minus_B",
        "hap_sets": "hapA1,hapA2 vs hapB1,hapB2",
        "genome_bp": group_rows[0]["genome_bp"] - group_rows[1]["genome_bp"],
        "gff_gene": a_gene - b_gene,
        "pep_records": group_rows[0]["pep_records"] - group_rows[1]["pep_records"],
        "cds_bp": group_rows[0]["cds_bp"] - group_rows[1]["cds_bp"],
    })
    write_tsv(out_dir / "AB_gene_number_difference.tsv", group_rows, [
        "group", "hap_sets", "genome_bp", "gff_gene", "pep_records", "cds_bp",
    ])

    busco = find_busco_summaries(busco_dir)
    busco_rows = []
    for hap in HAPS:
        row = {"hap_set": hap}
        row.update(busco.get(hap, {
            "busco_path": "",
            "complete_pct": "",
            "single_pct": "",
            "duplicated_pct": "",
            "fragmented_pct": "",
            "missing_pct": "",
            "n_buscos": "",
            "status": "missing",
        }))
        busco_rows.append(row)
    write_tsv(out_dir / "busco_summary.tsv", busco_rows, [
        "hap_set", "status", "complete_pct", "single_pct", "duplicated_pct",
        "fragmented_pct", "missing_pct", "n_buscos", "busco_path",
    ])

    write_tsv(out_dir / "synteny_summary.tsv", read_synteny_summaries(work_dir), [
        "comparison", "chromosome_pairs_with_alignment", "total_aligned_block_bp",
        "total_target_covered_bp", "total_query_covered_bp", "summary_path",
    ])

    rbh_summary = []
    loss_rows = []
    for query_hap, target_hap in PAIRWISE_RBH:
        rbh, query_loss, target_loss = rbh_for_pair(blast_dir, query_hap, target_hap, protein_ids)
        q_total = len(protein_ids.get(query_hap, []))
        t_total = len(protein_ids.get(target_hap, []))
        rbh_summary.append({
            "comparison": "%s_vs_%s" % (query_hap, target_hap),
            "query_hap": query_hap,
            "target_hap": target_hap,
            "query_total_proteins": q_total,
            "target_total_proteins": t_total,
            "reciprocal_best_hit_pairs": len(rbh),
            "query_without_RBH": len(query_loss),
            "target_without_RBH": len(target_loss),
            "query_retention_pct": "%.3f" % (100 * len(rbh) / q_total if q_total else 0),
            "target_retention_pct": "%.3f" % (100 * len(rbh) / t_total if t_total else 0),
        })
        for gene in query_loss:
            loss_rows.append({
                "comparison": "%s_vs_%s" % (query_hap, target_hap),
                "hap_set": query_hap,
                "gene_id": gene,
                "candidate_status": "no_reciprocal_best_hit_in_%s" % target_hap,
                "interpretation": "candidate retention/loss difference; not proof of biological gene loss",
            })
        for gene in target_loss:
            loss_rows.append({
                "comparison": "%s_vs_%s" % (query_hap, target_hap),
                "hap_set": target_hap,
                "gene_id": gene,
                "candidate_status": "no_reciprocal_best_hit_in_%s" % query_hap,
                "interpretation": "candidate retention/loss difference; not proof of biological gene loss",
            })
    write_tsv(out_dir / "rbh_gene_retention_summary.tsv", rbh_summary, [
        "comparison", "query_hap", "target_hap", "query_total_proteins", "target_total_proteins",
        "reciprocal_best_hit_pairs", "query_without_RBH", "target_without_RBH",
        "query_retention_pct", "target_retention_pct",
    ])
    write_tsv(out_dir / "gene_loss_candidates.tsv", loss_rows, [
        "comparison", "hap_set", "gene_id", "candidate_status", "interpretation",
    ])

    rank_rows = []
    for hap in HAPS:
        score, method = score_hap_for_representative(hap, metrics, busco)
        row = {
            "hap_set": hap,
            "group": "A" if hap in GROUPS["A"] else "B",
            "score": "%.4f" % score,
            "score_basis": method,
            "gene_count": metrics[hap]["gff_gene"],
            "genome_bp": metrics[hap]["genome_bp"],
            "busco_complete_pct": busco.get(hap, {}).get("complete_pct", ""),
            "busco_missing_pct": busco.get(hap, {}).get("missing_pct", ""),
        }
        rank_rows.append(row)
    chosen_a = sorted([r for r in rank_rows if r["group"] == "A"], key=lambda r: float(r["score"]), reverse=True)[0]["hap_set"]
    chosen_b = sorted([r for r in rank_rows if r["group"] == "B"], key=lambda r: float(r["score"]), reverse=True)[0]["hap_set"]
    for row in rank_rows:
        row["chosen_for_AplusB"] = "yes" if row["hap_set"] in {chosen_a, chosen_b} else "no"
    write_tsv(out_dir / "hap_representative_ranking.tsv", rank_rows, [
        "hap_set", "group", "score", "score_basis", "gene_count", "genome_bp",
        "busco_complete_pct", "busco_missing_pct", "chosen_for_AplusB",
    ])

    rec = [{
        "question": "Which haplotype should be used for species-level comparative genomics?",
        "answer": "%s+%s representative n=16 A+B set" % (chosen_a, chosen_b),
        "recommended_files": str(out_dir / "downstream_representative_AplusB"),
        "reason": "PZI has inferred A and B subgenome groups. A species-level representative should include one A haplotype and one B haplotype; hapA1 alone contains only the A group and is suitable only for A-subgenome-specific analyses.",
        "caveat": "A/B labels are inferred without parental data; retain FINAL_HAPLOTYPE_OVERVIEW.tsv and review notes with downstream results.",
    }]
    write_tsv(out_dir / "downstream_representative_recommendation.tsv", rec, [
        "question", "answer", "recommended_files", "reason", "caveat",
    ])
    write_representative_package(package_dir, out_dir, chosen_a, chosen_b)

    # The renamed package keeps its own manifest and chromosome map, while
    # the legacy phase-review tables remain in haploid/final_haplotypes.
    metadata_names = (
        "FINAL_HAPLOTYPE_OVERVIEW.tsv",
        "FINAL_REVIEW_NOTES.tsv",
        "chromosome_rename_map.tsv",
        "manifest.tsv",
        "validation.tsv",
    )
    legacy_metadata_dir = package_dir.parent / "haploid" / "final_haplotypes"
    for name in metadata_names:
        source = next(
            (
                candidate
                for candidate in (package_dir / name, legacy_metadata_dir / name)
                if candidate.is_file() and candidate.stat().st_size > 0
            ),
            None,
        )
        if source is not None:
            shutil.copy2(source, out_dir / name)

    (out_dir / "README.md").write_text(
        "Internal haplotype/subgenome comparison summary.\n\n"
        "Main outputs:\n"
        "- synteny_summary.tsv: minimap2 genome-synteny summary.\n"
        "- rbh_gene_retention_summary.tsv: reciprocal-best-hit protein retention summary.\n"
        "- gene_loss_candidates.tsv: genes without RBH in a counterpart haplotype; candidate list only.\n"
        "- AB_gene_number_difference.tsv: A-group vs B-group gene/genome size differences.\n"
        "- busco_summary.tsv: BUSCO protein-mode summaries if BUSCO was run.\n"
        "- downstream_representative_AplusB/: recommended A+B n=16 files for species-level comparative genomics.\n",
        encoding="utf-8",
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")

    paf = subparsers.add_parser("paf-summary")
    paf.add_argument("--paf", required=True, type=Path)
    paf.add_argument("--out", required=True, type=Path)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--package-dir", required=True, type=Path)
    collect_parser.add_argument("--work-dir", required=True, type=Path)
    collect_parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "paf-summary":
        write_paf_summary(args.paf, args.out)
    elif args.command == "collect":
        collect(args)
    else:
        raise SystemExit("ERROR: a subcommand is required: paf-summary or collect")


if __name__ == "__main__":
    main()
