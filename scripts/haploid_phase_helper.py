#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
from typing import Optional


PALETTE = [
    "#F6AD55", "#68D391", "#63B3ED", "#F687B3",
    "#4FD1C5", "#ED8936", "#9F7AEA", "#ECC94B",
    "#F56565", "#4299E1", "#48BB78", "#A0AEC0",
    "#B794F4", "#FC8181", "#76E4F7", "#FBD38D",
]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_fasta_headers(path: Path):
    headers = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip())
    return headers


def read_main_chromosomes_from_assembly_report(path: Path):
    rows = []
    if not path.is_file():
        return rows
    header = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("# Sequence-Name\t"):
                header = line.lstrip("# ").split("\t")
                continue
            if not header or line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            if (
                row.get("Sequence-Role") == "assembled-molecule"
                and row.get("Assigned-Molecule-Location/Type") in {"Chromosome", "Linkage Group"}
            ):
                rows.append({
                    "sequence_name": row.get("Sequence-Name", ""),
                    "assigned_molecule": row.get("Assigned-Molecule", ""),
                    "genbank_accn": row.get("GenBank-Accn", ""),
                    "refseq_accn": row.get("RefSeq-Accn", ""),
                    "length_bp": row.get("Sequence-Length", ""),
                    "assembly_unit": row.get("Assembly-Unit", ""),
                })
    return rows


def main_reference_order(fasta_headers, report_path: Path, expected_count: Optional[int] = None):
    fasta_names = [h.split()[0] for h in fasta_headers]
    fasta_set = set(fasta_names)
    report_rows = read_main_chromosomes_from_assembly_report(report_path)
    order = []
    manifest = []
    for row in report_rows:
        candidates = [row["genbank_accn"], row["refseq_accn"], row["sequence_name"]]
        chosen = next((name for name in candidates if name and name != "na" and name in fasta_set), "")
        if not chosen:
            continue
        order.append(chosen)
        manifest.append({
            "fasta_chr": chosen,
            "sequence_name": row["sequence_name"],
            "assigned_molecule": row["assigned_molecule"],
            "genbank_accn": row["genbank_accn"],
            "refseq_accn": row["refseq_accn"],
            "length_bp": row["length_bp"],
            "assembly_unit": row["assembly_unit"],
        })
    if order:
        if expected_count and len(order) > expected_count:
            order = order[:expected_count]
            manifest = manifest[:expected_count]
        return order, manifest
    fallback = [
        name for name in fasta_names
        if not any(token in name.lower() for token in ["scaffold", "contig", "unplaced"])
    ]
    if not fallback:
        fallback = fasta_names
    manifest = [{
        "fasta_chr": name,
        "sequence_name": name,
        "assigned_molecule": "",
        "genbank_accn": "",
        "refseq_accn": "",
        "length_bp": "",
        "assembly_unit": "fallback_from_fasta_order",
    } for name in fallback]
    if expected_count and len(fallback) > expected_count:
        fallback = fallback[:expected_count]
        manifest = manifest[:expected_count]
    return fallback, manifest


def read_tsv(path: Path):
    if not path.is_file():
        raise SystemExit(f"ERROR: missing file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def natural_key(text: str):
    parts = []
    buf = ""
    is_digit = False
    for ch in text:
        if ch.isdigit() != is_digit and buf:
            parts.append(int(buf) if is_digit else buf)
            buf = ch
            is_digit = ch.isdigit()
        else:
            buf += ch
            is_digit = ch.isdigit()
    if buf:
        parts.append(int(buf) if is_digit else buf)
    return parts


def parse_paf(path: Path, source_label: str, min_mapq: int, min_aln_len: int):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])
            if mapq < min_mapq or aln_len < min_aln_len:
                continue
            rows.append({
                "source": source_label,
                "query": qname,
                "query_len": qlen,
                "query_start": qstart,
                "query_end": qend,
                "strand": strand,
                "target": tname,
                "target_len": tlen,
                "target_start": tstart,
                "target_end": tend,
                "matches": matches,
                "aln_len": aln_len,
                "mapq": mapq,
                "identity": matches / aln_len if aln_len else 0.0,
            })
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


def summarize_pzi_vs_ref(rows, pzi_order, ref_order, ref_label):
    by_pzi = defaultdict(lambda: defaultdict(lambda: {
        "covered": 0,
        "matches": 0,
        "aln_len": 0,
        "plus": 0,
        "minus": 0,
        "intervals": [],
    }))
    pzi_lengths = {}
    ref_lengths = {}
    ref_col = f"{ref_label}_chr"
    ref_len_col = f"{ref_label}_length_bp"

    for row in rows:
        pzi = row["query"]
        ref = row["target"]
        pzi_lengths[pzi] = row["query_len"]
        ref_lengths[ref] = row["target_len"]
        by_pzi[pzi][ref]["matches"] += row["matches"]
        by_pzi[pzi][ref]["aln_len"] += row["aln_len"]
        by_pzi[pzi][ref]["intervals"].append((row["query_start"], row["query_end"]))
        if row["strand"] == "+":
            by_pzi[pzi][ref]["plus"] += row["aln_len"]
        else:
            by_pzi[pzi][ref]["minus"] += row["aln_len"]

    out = []
    best_rows = []
    for pzi in sorted(pzi_order, key=natural_key):
        best_ref = ""
        best_cov = -1
        best_row = None
        for ref, stats in by_pzi.get(pzi, {}).items():
            cov = interval_union_bp(stats["intervals"])
            if cov > best_cov:
                best_cov = cov
                best_ref = ref
        total_cov = 0
        total_aln = 0
        best_row = None
        for ref in ref_order:
            stats = by_pzi.get(pzi, {}).get(ref)
            if not stats:
                continue
            covered_bp = interval_union_bp(stats["intervals"])
            aln_len = stats["aln_len"]
            identity = stats["matches"] / aln_len if aln_len else 0.0
            orientation = "+" if stats["plus"] >= stats["minus"] else "-"
            out.append({
                "pzi_chr": pzi,
                ref_col: ref,
                "pzi_length_bp": str(pzi_lengths.get(pzi, "")),
                ref_len_col: str(ref_lengths.get(ref, "")),
                "aligned_block_bp": str(aln_len),
                "covered_bp": str(covered_bp),
                "coverage_ratio": f"{(covered_bp / pzi_lengths.get(pzi, 1)):.6f}" if pzi_lengths.get(pzi) else "",
                "identity": f"{identity:.6f}",
                "orientation": orientation,
                "best_hit_for_pzi_chr": "yes" if ref == best_ref else "no",
            })
            total_cov += covered_bp
            total_aln += aln_len
            if ref == best_ref:
                best_row = out[-1]
        if best_row is None:
            best_row = {
                "pzi_chr": pzi,
                ref_col: "",
                "pzi_length_bp": str(pzi_lengths.get(pzi, "")),
                ref_len_col: "",
                "aligned_block_bp": "0",
                "covered_bp": "0",
                "coverage_ratio": "0.000000",
                "identity": "0.000000",
                "orientation": "",
                "best_hit_for_pzi_chr": "yes",
            }
        best_rows.append(best_row)
    return out, best_rows


def summarize_pzi_vs_pcam(rows, pzi_order, pcam_order):
    return summarize_pzi_vs_ref(rows, pzi_order, pcam_order, "pcam")


def summarize_pzi_vs_ppa(rows, pzi_order, ppa_order):
    return summarize_pzi_vs_ref(rows, pzi_order, ppa_order, "ppa")


def parse_reference_spec(text: str):
    parts = text.split("|")
    if len(parts) not in {4, 5}:
        raise SystemExit(
            "ERROR: --n8-reference must use label|fasta|assembly_report|paf or "
            "label|fasta|assembly_report|paf|expected_count format, "
            f"observed: {text}"
        )
    label, fasta, report, paf = parts[:4]
    if not label:
        raise SystemExit(f"ERROR: empty reference label in --n8-reference {text}")
    expected_count = 8
    if len(parts) == 5:
        try:
            expected_count = int(parts[4])
        except ValueError:
            raise SystemExit(f"ERROR: expected_count must be an integer in --n8-reference {text}")
    return {
        "label": label,
        "fasta": Path(fasta),
        "assembly_report": Path(report),
        "paf": Path(paf),
        "expected_count": expected_count,
    }


def parse_chr_list(text: str):
    return [item.strip() for item in text.split(",") if item.strip()]


def build_homolog_quartets(pzi_rows, provenance_rows, quartet_map):
    by_pcam = defaultdict(list)
    for row in pzi_rows:
        if row["best_hit_for_pzi_chr"] == "yes" and row["pcam_chr"]:
            by_pcam[row["pcam_chr"]].append(row)

    copy_lookup = {row["chr"]: row for row in provenance_rows}

    quartets = []
    for idx, pcam_chr in enumerate(sorted(by_pcam, key=natural_key), start=1):
        members = by_pcam[pcam_chr]
        members.sort(key=lambda r: (-float(r["covered_bp"] or 0), -float(r["identity"] or 0)))
        while len(members) < 4:
            members.append({
                "pzi_chr": "",
                "pcam_chr": pcam_chr,
                "aligned_block_bp": "0",
                "covered_bp": "0",
                "identity": "0",
                "orientation": "",
                "best_hit_for_pzi_chr": "no",
            })
        quartet_id = f"quartet{idx:02d}"
        quartets.append({
            "quartet_id": quartet_id,
            "pcam_chr": pcam_chr,
            "pzi_chr_1": members[0]["pzi_chr"],
            "pzi_chr_2": members[1]["pzi_chr"],
            "pzi_chr_3": members[2]["pzi_chr"],
            "pzi_chr_4": members[3]["pzi_chr"],
        })
        quartet_map[quartet_id] = [m["pzi_chr"] for m in members if m["pzi_chr"]]
    quartet_rows = []
    for quartet in quartets:
        quartet_id = quartet["quartet_id"]
        for order, chr_name in enumerate([quartet["pzi_chr_1"], quartet["pzi_chr_2"], quartet["pzi_chr_3"], quartet["pzi_chr_4"]], start=1):
            if not chr_name:
                continue
            quartet_rows.append({
                "quartet_id": quartet_id,
                "quartet_rank": str(order),
                "pcam_chr": quartet["pcam_chr"],
                "pzi_chr": chr_name,
                "source_copy": copy_lookup.get(chr_name, {}).get("source_copy", ""),
                "dominant_hap": copy_lookup.get(chr_name, {}).get("dominant_hap", ""),
                "copy_provenance_status": copy_lookup.get(chr_name, {}).get("status", ""),
                "copy_dominance": copy_lookup.get(chr_name, {}).get("dominance", ""),
                "copy_top_query": copy_lookup.get(chr_name, {}).get("top_hap1_query", "") if copy_lookup.get(chr_name, {}).get("dominant_hap", "") == "hap1" else copy_lookup.get(chr_name, {}).get("top_hap2_query", ""),
            })
    return quartets, quartet_rows


def best_index(best_rows, ref_label):
    ref_col = f"{ref_label}_chr"
    len_col = f"{ref_label}_length_bp"
    return {
        row["pzi_chr"]: {
            "ref_chr": row.get(ref_col, ""),
            "ref_length_bp": row.get(len_col, ""),
            "covered_bp": row.get("covered_bp", "0"),
            "aligned_block_bp": row.get("aligned_block_bp", "0"),
            "coverage_ratio": row.get("coverage_ratio", "0"),
            "identity": row.get("identity", "0"),
            "orientation": row.get("orientation", ""),
        }
        for row in best_rows
    }


def build_reference_review(best_rows, ref_label):
    ref_col = f"{ref_label}_chr"
    len_col = f"{ref_label}_length_bp"
    grouped = defaultdict(list)
    for row in best_rows:
        if row.get(ref_col):
            grouped[row[ref_col]].append(row)
    out = []
    for ref_chr in sorted(grouped, key=natural_key):
        members = grouped[ref_chr]
        pzi_total = sum(int(float(row.get("pzi_length_bp", "0") or 0)) for row in members)
        ref_len = 0
        for row in members:
            if row.get(len_col):
                ref_len = int(float(row[len_col]))
                break
        ratio = pzi_total / ref_len if ref_len else 0.0
        count = len(members)
        if count == 4:
            event = "balanced_quartet"
            eej_ncf_rct = "none"
            status = "pass"
            note = "Four PZI chromosomes have this reference chromosome as best hit."
        elif count >= 5:
            event = "possible_fusion_or_complex_rearrangement"
            eej_ncf_rct = "NCF_like"
            status = "review"
            note = "More than four PZI chromosomes point here; inspect for nested chromosome fusion or fragmented assignment."
        elif count == 3:
            event = "possible_break_or_missing_partner"
            eej_ncf_rct = "EEJ_like"
            status = "review"
            note = "Three PZI chromosomes point here; inspect for breakpoint, missing partner, or weak synteny."
        elif count == 2:
            event = "possible_fusion_or_split"
            eej_ncf_rct = "RCT_like"
            status = "review"
            note = "Two PZI chromosomes point here; inspect for fusion/split-like structure."
        else:
            event = "weak_or_singleton_assignment"
            eej_ncf_rct = "complex"
            status = "review"
            note = "Singleton or weak assignment; check dotplot and heatmap."
        if ref_len and ratio >= 4.8:
            event += "; length_excess"
            note += " Combined PZI length is very large relative to the reference chromosome."
        elif ref_len and 0 < ratio <= 1.2 and count <= 3:
            event += "; length_short"
            note += " Combined PZI length is short relative to the reference chromosome."
        out.append({
            "reference_species": ref_label,
            "reference_chr": ref_chr,
            "pzi_members": ",".join(sorted([row["pzi_chr"] for row in members], key=natural_key)),
            "member_count": str(count),
            "reference_length_bp": str(ref_len),
            "pzi_member_total_bp": str(pzi_total),
            "length_ratio_to_reference": f"{ratio:.4f}",
            "event_class": event,
            "eej_ncf_rct_candidate": eej_ncf_rct,
            "review_status": status,
            "review_note": note,
        })
    return out


def build_trispecies_pzi_table(pzi_order, pcam_best, ppa_best):
    pcam = best_index(pcam_best, "pcam")
    ppa = best_index(ppa_best, "ppa")
    rows = []
    for pzi_chr in sorted(pzi_order, key=natural_key):
        pcam_row = pcam.get(pzi_chr, {})
        ppa_row = ppa.get(pzi_chr, {})
        pcam_cov = float(pcam_row.get("coverage_ratio", "0") or 0)
        ppa_cov = float(ppa_row.get("coverage_ratio", "0") or 0)
        if pcam_cov < 0.05 or ppa_cov < 0.05:
            status = "review"
            note = "At least one reference gives weak best-hit coverage."
        elif pcam_row.get("orientation") and ppa_row.get("orientation") and pcam_row.get("orientation") != ppa_row.get("orientation"):
            status = "review"
            note = "PCAM and PPA best-hit orientations differ; inspect for inversion or rearranged segment."
        else:
            status = "pass"
            note = "PCAM and PPA best-hit evidence is broadly consistent."
        rows.append({
            "pzi_chr": pzi_chr,
            "pcam_chr": pcam_row.get("ref_chr", ""),
            "pcam_coverage_ratio": pcam_row.get("coverage_ratio", "0"),
            "pcam_orientation": pcam_row.get("orientation", ""),
            "ppa_chr": ppa_row.get("ref_chr", ""),
            "ppa_coverage_ratio": ppa_row.get("coverage_ratio", "0"),
            "ppa_orientation": ppa_row.get("orientation", ""),
            "review_status": status,
            "review_note": note,
        })
    return rows


def group_pzi_by_ref(best_rows, ref_label):
    ref_col = f"{ref_label}_chr"
    grouped = defaultdict(list)
    for row in best_rows:
        ref_chr = row.get(ref_col, "")
        pzi_chr = row.get("pzi_chr", "")
        if ref_chr and pzi_chr:
            grouped[ref_chr].append(pzi_chr)
    return {
        ref_chr: sorted(set(chrs), key=natural_key)
        for ref_chr, chrs in grouped.items()
    }


def build_n8_reference_outputs(args, out_dir, specs, pzi_chr_headers, chr_order, min_mapq, min_aln_len):
    outputs = {}
    manifest_rows = []
    for spec in specs:
        label = spec["label"]
        ref_headers, ref_manifest = main_reference_order(
            read_fasta_headers(spec["fasta"]),
            spec["assembly_report"],
            spec.get("expected_count", 8),
        )
        for row in ref_manifest:
            row = dict(row)
            row["reference_label"] = label
            manifest_rows.append(row)
        paf_rows = parse_paf(spec["paf"], label, min_mapq, min_aln_len)
        summary, best = summarize_pzi_vs_ref(paf_rows, pzi_chr_headers, ref_headers, label)
        write_tsv(out_dir / f"pzi_vs_{label}_chr_alignment.tsv", summary, [
            "pzi_chr", f"{label}_chr", "pzi_length_bp", f"{label}_length_bp", "aligned_block_bp",
            "covered_bp", "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
        ])
        write_tsv(out_dir / f"pzi_vs_{label}_best_hits.tsv", best, [
            "pzi_chr", f"{label}_chr", "pzi_length_bp", f"{label}_length_bp", "aligned_block_bp",
            "covered_bp", "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
        ])
        write_tsv(out_dir / f"{label}_rearrangement_review.tsv", build_reference_review(best, label), [
            "reference_species", "reference_chr", "pzi_members", "member_count", "reference_length_bp",
            "pzi_member_total_bp", "length_ratio_to_reference", "event_class", "eej_ncf_rct_candidate",
            "review_status", "review_note"
        ])
        build_reference_query_summary(paf_rows, out_dir / f"pzi_vs_{label}_alignment_summary.tsv")
        outputs[label] = {
            "label": label,
            "headers": ref_headers,
            "manifest": ref_manifest,
            "paf_rows": paf_rows,
            "summary": summary,
            "best": best,
        }
    if manifest_rows:
        write_tsv(out_dir / "n8_reference_main_chromosomes.tsv", manifest_rows, [
            "reference_label", "fasta_chr", "sequence_name", "assigned_molecule", "genbank_accn",
            "refseq_accn", "length_bp", "assembly_unit"
        ])
    return outputs


def build_outgroup_outputs(out_dir, specs, pzi_chr_headers, focus_chrs, min_mapq, min_aln_len, bin_size):
    outputs = {}
    manifest_rows = []
    segment_rows = []
    chr_vote_rows = []
    bin_vote_rows = []
    focus_set = set(focus_chrs)

    for spec in specs:
        label = spec["label"]
        ref_headers, ref_manifest = main_reference_order(
            read_fasta_headers(spec["fasta"]),
            spec["assembly_report"],
            spec.get("expected_count", None),
        )
        for row in ref_manifest:
            row = dict(row)
            row["reference_label"] = label
            manifest_rows.append(row)

        paf_rows = parse_paf(spec["paf"], label, min_mapq, min_aln_len)
        summary, best = summarize_pzi_vs_ref(paf_rows, pzi_chr_headers, ref_headers, label)
        write_tsv(out_dir / f"pzi_vs_{label}_chr_alignment.tsv", summary, [
            "pzi_chr", f"{label}_chr", "pzi_length_bp", f"{label}_length_bp", "aligned_block_bp",
            "covered_bp", "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
        ])
        write_tsv(out_dir / f"pzi_vs_{label}_best_hits.tsv", best, [
            "pzi_chr", f"{label}_chr", "pzi_length_bp", f"{label}_length_bp", "aligned_block_bp",
            "covered_bp", "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
        ])
        write_tsv(out_dir / f"{label}_rearrangement_review.tsv", build_reference_review(best, label), [
            "reference_species", "reference_chr", "pzi_members", "member_count", "reference_length_bp",
            "pzi_member_total_bp", "length_ratio_to_reference", "event_class", "eej_ncf_rct_candidate",
            "review_status", "review_note"
        ])
        build_reference_query_summary(paf_rows, out_dir / f"pzi_vs_{label}_alignment_summary.tsv")

        chr_votes = defaultdict(lambda: {
            "aligned_bp": 0,
            "matches": 0,
            "segments": 0,
            "plus_bp": 0,
            "minus_bp": 0,
            "query_len": 0,
            "target_len": 0,
        })
        bin_votes = defaultdict(lambda: {
            "aligned_bp": 0,
            "matches": 0,
            "segments": 0,
            "plus_bp": 0,
            "minus_bp": 0,
            "query_len": 0,
            "target_len": 0,
        })

        for row in paf_rows:
            if focus_set and row["query"] not in focus_set:
                continue
            segment_rows.append({
                "outgroup": label,
                "pzi_chr": row["query"],
                "pzi_start": str(row["query_start"]),
                "pzi_end": str(row["query_end"]),
                "outgroup_chr": row["target"],
                "outgroup_start": str(row["target_start"]),
                "outgroup_end": str(row["target_end"]),
                "strand": row["strand"],
                "aligned_bp": str(row["aln_len"]),
                "matching_bp": str(row["matches"]),
                "identity": f"{row['identity']:.6f}",
                "mapq": str(row["mapq"]),
            })
            chr_key = (label, row["query"], row["target"])
            chr_votes[chr_key]["aligned_bp"] += row["aln_len"]
            chr_votes[chr_key]["matches"] += row["matches"]
            chr_votes[chr_key]["segments"] += 1
            chr_votes[chr_key]["query_len"] = row["query_len"]
            chr_votes[chr_key]["target_len"] = row["target_len"]
            if row["strand"] == "+":
                chr_votes[chr_key]["plus_bp"] += row["aln_len"]
            else:
                chr_votes[chr_key]["minus_bp"] += row["aln_len"]

            start_bin = row["query_start"] // bin_size
            end_bin = max(start_bin, (max(row["query_end"] - 1, row["query_start"])) // bin_size)
            for bin_id in range(start_bin, end_bin + 1):
                bin_start = bin_id * bin_size
                bin_end = bin_start + bin_size
                overlap = max(0, min(row["query_end"], bin_end) - max(row["query_start"], bin_start))
                if overlap <= 0:
                    continue
                bin_key = (label, row["query"], str(bin_id), row["target"])
                bin_votes[bin_key]["aligned_bp"] += overlap
                bin_votes[bin_key]["matches"] += int(row["matches"] * (overlap / row["aln_len"])) if row["aln_len"] else 0
                bin_votes[bin_key]["segments"] += 1
                bin_votes[bin_key]["query_len"] = row["query_len"]
                bin_votes[bin_key]["target_len"] = row["target_len"]
                if row["strand"] == "+":
                    bin_votes[bin_key]["plus_bp"] += overlap
                else:
                    bin_votes[bin_key]["minus_bp"] += overlap

        total_by_chr = defaultdict(int)
        for label2, pzi_chr, target_chr in chr_votes:
            total_by_chr[(label2, pzi_chr)] += chr_votes[(label2, pzi_chr, target_chr)]["aligned_bp"]
        for (label2, pzi_chr, target_chr), stats in sorted(chr_votes.items(), key=lambda kv: (natural_key(kv[0][0]), natural_key(kv[0][1]), -kv[1]["aligned_bp"])):
            total = total_by_chr[(label2, pzi_chr)]
            chr_vote_rows.append({
                "outgroup": label2,
                "pzi_chr": pzi_chr,
                "outgroup_chr": target_chr,
                "aligned_bp": str(stats["aligned_bp"]),
                "vote_fraction": f"{(stats['aligned_bp'] / total):.6f}" if total else "0.000000",
                "segments": str(stats["segments"]),
                "orientation": "+" if stats["plus_bp"] >= stats["minus_bp"] else "-",
                "identity": f"{(stats['matches'] / stats['aligned_bp']):.6f}" if stats["aligned_bp"] else "0.000000",
                "pzi_length_bp": str(stats["query_len"]),
                "outgroup_length_bp": str(stats["target_len"]),
            })

        total_by_bin = defaultdict(int)
        for label2, pzi_chr, bin_id, target_chr in bin_votes:
            total_by_bin[(label2, pzi_chr, bin_id)] += bin_votes[(label2, pzi_chr, bin_id, target_chr)]["aligned_bp"]
        for (label2, pzi_chr, bin_id, target_chr), stats in sorted(bin_votes.items(), key=lambda kv: (natural_key(kv[0][0]), natural_key(kv[0][1]), int(kv[0][2]), -kv[1]["aligned_bp"])):
            total = total_by_bin[(label2, pzi_chr, bin_id)]
            bin_start = int(bin_id) * bin_size
            bin_vote_rows.append({
                "outgroup": label2,
                "pzi_chr": pzi_chr,
                "bin_id": bin_id,
                "bin_start": str(bin_start),
                "bin_end": str(bin_start + bin_size),
                "outgroup_chr": target_chr,
                "aligned_bp": str(stats["aligned_bp"]),
                "vote_fraction": f"{(stats['aligned_bp'] / total):.6f}" if total else "0.000000",
                "segments": str(stats["segments"]),
                "orientation": "+" if stats["plus_bp"] >= stats["minus_bp"] else "-",
                "identity": f"{(stats['matches'] / stats['aligned_bp']):.6f}" if stats["aligned_bp"] else "0.000000",
            })

        outputs[label] = {
            "label": label,
            "headers": ref_headers,
            "manifest": ref_manifest,
            "paf_rows": paf_rows,
            "summary": summary,
            "best": best,
        }

    if manifest_rows:
        write_tsv(out_dir / "outgroup_main_chromosomes.tsv", manifest_rows, [
            "reference_label", "fasta_chr", "sequence_name", "assigned_molecule", "genbank_accn",
            "refseq_accn", "length_bp", "assembly_unit"
        ])
    write_tsv(out_dir / "outgroup_segments.tsv", segment_rows, [
        "outgroup", "pzi_chr", "pzi_start", "pzi_end", "outgroup_chr", "outgroup_start",
        "outgroup_end", "strand", "aligned_bp", "matching_bp", "identity", "mapq"
    ])
    write_tsv(out_dir / "outgroup_chr_votes.tsv", chr_vote_rows, [
        "outgroup", "pzi_chr", "outgroup_chr", "aligned_bp", "vote_fraction", "segments",
        "orientation", "identity", "pzi_length_bp", "outgroup_length_bp"
    ])
    write_tsv(out_dir / "outgroup_bin_votes.tsv", bin_vote_rows, [
        "outgroup", "pzi_chr", "bin_id", "bin_start", "bin_end", "outgroup_chr",
        "aligned_bp", "vote_fraction", "segments", "orientation", "identity"
    ])
    return outputs


def build_multi_reference_consensus(chr_order, ppa_best, pcam_best, n8_outputs):
    ppa_groups = group_pzi_by_ref(ppa_best, "ppa")
    pcam_idx = best_index(pcam_best, "pcam")
    n8_idx = {label: best_index(data["best"], label) for label, data in n8_outputs.items()}
    ppa_ref_col = "ppa_chr"
    rows = []
    for ppa_chr in sorted(ppa_groups, key=natural_key):
        members = ppa_groups[ppa_chr]
        pcam_votes = defaultdict(int)
        n8_vote_parts = []
        consensus_votes = defaultdict(int)
        for member in members:
            pcam_chr = pcam_idx.get(member, {}).get("ref_chr", "")
            if pcam_chr:
                pcam_votes[pcam_chr] += 1
            for label, idx in n8_idx.items():
                ref_chr = idx.get(member, {}).get("ref_chr", "")
                cov = idx.get(member, {}).get("coverage_ratio", "0")
                if ref_chr:
                    consensus_votes[f"{label}:{ref_chr}"] += 1
                n8_vote_parts.append(f"{member}:{label}:{ref_chr}:{cov}")
        top_votes = sorted(consensus_votes.items(), key=lambda kv: (-kv[1], natural_key(kv[0])))
        top_n8 = top_votes[0][0] if top_votes else ""
        top_n8_count = top_votes[0][1] if top_votes else 0
        dominant_pcam = sorted(pcam_votes.items(), key=lambda kv: (-kv[1], natural_key(kv[0])))[0][0] if pcam_votes else ""
        status = "pass" if top_n8 and top_n8_count >= max(2, len(members)) and len(set(members)) == 2 else "review"
        rows.append({
            "ppa_chr": ppa_chr,
            "pzi_members": ",".join(members),
            "member_count": str(len(members)),
            "pcam_vote": dominant_pcam,
            "n8_consensus_vote": top_n8,
            "n8_consensus_count": str(top_n8_count),
            "all_n8_votes": ";".join(n8_vote_parts),
            "status": status,
            "review_note": "Use n8 consensus to pair PPA-derived two-chromosome groups into four-chromosome homologous quartets.",
        })
    return rows


def build_consensus_quartets(consensus_rows):
    grouped = defaultdict(list)
    for row in consensus_rows:
        key = row["n8_consensus_vote"] or row["pcam_vote"] or row["ppa_chr"]
        grouped[key].append(row)
    quartet_rows = []
    member_rows = []
    for idx, key in enumerate(sorted(grouped, key=natural_key), start=1):
        rows = grouped[key]
        members = []
        ppa_groups = []
        for row in rows:
            ppa_groups.append(row["ppa_chr"])
            members.extend([x for x in row["pzi_members"].split(",") if x])
        members = sorted(set(members), key=natural_key)
        quartet_id = f"consensus_quartet{idx:02d}"
        padded = members + [""] * max(0, 4 - len(members))
        status = "pass" if len(members) == 4 and len(rows) == 2 else "review"
        quartet_rows.append({
            "consensus_quartet_id": quartet_id,
            "consensus_key": key,
            "ppa_groups": ",".join(sorted(ppa_groups, key=natural_key)),
            "member_count": str(len(members)),
            "pzi_chr_1": padded[0],
            "pzi_chr_2": padded[1],
            "pzi_chr_3": padded[2],
            "pzi_chr_4": padded[3],
            "status": status,
            "review_note": "Expected pass: two PPA groups combine into four PZI chromosomes under the same n8 consensus.",
        })
        for rank, chr_name in enumerate(members, start=1):
            member_rows.append({
                "consensus_quartet_id": quartet_id,
                "consensus_key": key,
                "rank": str(rank),
                "pzi_chr": chr_name,
                "ppa_groups": ",".join(sorted(ppa_groups, key=natural_key)),
                "status": status,
            })
    return quartet_rows, member_rows


def build_color_maps_from_quartets(quartets, quartet_map):
    color_by_pzi = {}
    color_by_pcam = {}
    for idx, quartet in enumerate(sorted(quartets, key=lambda r: natural_key(r["quartet_id"]))):
        color = PALETTE[idx % len(PALETTE)]
        color_by_pcam[quartet["pcam_chr"]] = color
        for pzi_chr in quartet_map.get(quartet["quartet_id"], []):
            color_by_pzi[pzi_chr] = color
    return color_by_pzi, color_by_pcam


def build_reference_color_map(best_rows, ref_label, color_by_pzi):
    ref_col = f"{ref_label}_chr"
    votes = defaultdict(lambda: defaultdict(float))
    for row in best_rows:
        ref_chr = row.get(ref_col, "")
        pzi_chr = row.get("pzi_chr", "")
        if not ref_chr or not pzi_chr:
            continue
        color = color_by_pzi.get(pzi_chr)
        if color:
            votes[ref_chr][color] += float(row.get("covered_bp", "0") or 0)
    out = {}
    for ref_chr, color_votes in votes.items():
        out[ref_chr] = max(color_votes.items(), key=lambda kv: kv[1])[0]
    return out


def color_by_pzi_from_consensus(consensus_quartets):
    color_by_pzi = {}
    for idx, quartet in enumerate(sorted(consensus_quartets, key=lambda r: natural_key(r["consensus_quartet_id"]))):
        color = PALETTE[idx % len(PALETTE)]
        for key in ["pzi_chr_1", "pzi_chr_2", "pzi_chr_3", "pzi_chr_4"]:
            chr_name = quartet.get(key, "")
            if chr_name:
                color_by_pzi[chr_name] = color
    return color_by_pzi


def write_axis_colors(path: Path, color_by_pzi, color_by_ref, ref_label):
    rows = []
    for ref_chr, color in sorted(color_by_ref.items(), key=lambda kv: natural_key(kv[0])):
        rows.append({"axis": "target", "chr": ref_chr, "color": color, "label": f"{ref_label}_group"})
    for pzi_chr, color in sorted(color_by_pzi.items(), key=lambda kv: natural_key(kv[0])):
        rows.append({"axis": "query", "chr": pzi_chr, "color": color, "label": "pzi_quartet_group"})
    write_tsv(path, rows, ["axis", "chr", "color", "label"])


def write_axis_segments(path: Path, paf_rows, color_by_pzi, color_by_ref, min_segment_bp=50000):
    rows = []
    seen = set()
    for row in paf_rows:
        ref_chr = row["target"]
        pzi_chr = row["query"]
        if row["aln_len"] < min_segment_bp:
            continue
        color = color_by_pzi.get(pzi_chr) or color_by_ref.get(ref_chr) or "#CBD5E0"
        label = f"{pzi_chr}<->{ref_chr}"
        entries = [
            ("target", ref_chr, str(row["target_start"]), str(row["target_end"]), color, label),
            ("query", pzi_chr, str(row["query_start"]), str(row["query_end"]), color, label),
        ]
        for entry in entries:
            if entry in seen:
                continue
            seen.add(entry)
            rows.append({
                "axis": entry[0],
                "chr": entry[1],
                "start": entry[2],
                "end": entry[3],
                "color": entry[4],
                "label": entry[5],
            })
    write_tsv(path, rows, ["axis", "chr", "start", "end", "color", "label"])


def build_all_vs_all_matrix(rows, chr_order):
    matrix = defaultdict(int)
    for row in rows:
        if row["query"] in chr_order and row["target"] in chr_order:
            matrix[(row["query"], row["target"])] += row["aln_len"]
    out = []
    for q in chr_order:
        for t in chr_order:
            out.append({
                "query_chr": q,
                "target_chr": t,
                "aligned_bp": str(matrix.get((q, t), 0)),
            })
    return out


def infer_haplotype_candidates(quartet_map, provenance_rows, corrected_homolog_rows):
    prov = {row["chr"]: row for row in provenance_rows if row["chr"].startswith("chr")}
    pair_rows = []
    for row in corrected_homolog_rows:
        copy1_chr = row.get("copy1_chr", "") or row.get("copy1_source", "")
        copy2_chr = row.get("copy2_chr", "") or row.get("copy2_source", "")
        if copy1_chr and copy2_chr:
            pair_rows.append({
                "pair_id": row.get("pair_id", ""),
                "copy1_chr": copy1_chr,
                "copy2_chr": copy2_chr,
            })
    rows = []
    for quartet_id in sorted(quartet_map, key=natural_key):
        members = quartet_map[quartet_id]
        member_haps = [prov.get(chr_name, {}).get("dominant_hap", "") for chr_name in members]
        quartet_member_set = set(members)
        quartet_pairs = [
            row for row in pair_rows
            if row["copy1_chr"] in quartet_member_set and row["copy2_chr"] in quartet_member_set
        ]
        quartet_pairs.sort(key=lambda r: natural_key(r["pair_id"] or r["copy1_chr"]))

        if len(quartet_pairs) >= 2:
            a1 = quartet_pairs[0]["copy1_chr"]
            a2 = quartet_pairs[0]["copy2_chr"]
            b1 = quartet_pairs[1]["copy1_chr"]
            b2 = quartet_pairs[1]["copy2_chr"]
            status = "provisional_pair_based"
            candidate_basis = "two copy1-copy2 homolog pairs inside quartet"
            subgenome_label_warning = "A/B pair labels are provisional; confirm with dotplot/heatmap or subgenome-informative reference."
        else:
            copy1_members = [chr_name for chr_name in members if prov.get(chr_name, {}).get("source_copy", "") == "copy1"]
            copy2_members = [chr_name for chr_name in members if prov.get(chr_name, {}).get("source_copy", "") == "copy2"]
            a1 = copy1_members[0] if len(copy1_members) > 0 else (members[0] if len(members) > 0 else "")
            b1 = copy1_members[1] if len(copy1_members) > 1 else (members[1] if len(members) > 1 else "")
            a2 = copy2_members[0] if len(copy2_members) > 0 else (members[2] if len(members) > 2 else "")
            b2 = copy2_members[1] if len(copy2_members) > 1 else (members[3] if len(members) > 3 else "")
            status = "review"
            candidate_basis = "fallback copy split; homolog-pair evidence incomplete"
            subgenome_label_warning = "Review required; A/B pair labels are not resolved."
        row = {
            "quartet_id": quartet_id,
            "member_1": members[0] if len(members) > 0 else "",
            "member_2": members[1] if len(members) > 1 else "",
            "member_3": members[2] if len(members) > 2 else "",
            "member_4": members[3] if len(members) > 3 else "",
            "member_1_hap": member_haps[0] if len(member_haps) > 0 else "",
            "member_2_hap": member_haps[1] if len(member_haps) > 1 else "",
            "member_3_hap": member_haps[2] if len(member_haps) > 2 else "",
            "member_4_hap": member_haps[3] if len(member_haps) > 3 else "",
            "member_1_copy": prov.get(members[0], {}).get("source_copy", "") if len(members) > 0 else "",
            "member_2_copy": prov.get(members[1], {}).get("source_copy", "") if len(members) > 1 else "",
            "member_3_copy": prov.get(members[2], {}).get("source_copy", "") if len(members) > 2 else "",
            "member_4_copy": prov.get(members[3], {}).get("source_copy", "") if len(members) > 3 else "",
            "candidate_A1": a1,
            "candidate_A2": a2,
            "candidate_B1": b1,
            "candidate_B2": b2,
            "candidate_A1_copy": prov.get(a1, {}).get("source_copy", ""),
            "candidate_A2_copy": prov.get(a2, {}).get("source_copy", ""),
            "candidate_B1_copy": prov.get(b1, {}).get("source_copy", ""),
            "candidate_B2_copy": prov.get(b2, {}).get("source_copy", ""),
            "candidate_basis": candidate_basis,
            "subgenome_label_warning": subgenome_label_warning,
            "status": status,
        }
        rows.append(row)
    return rows


def build_source_copy_lookup(copy_map_rows):
    lookup = {}
    for row in copy_map_rows:
        copy1_chr = row.get("copy1_chr", "") or row.get("copy1_source", "")
        copy2_chr = row.get("copy2_chr", "") or row.get("copy2_source", "")
        if copy1_chr:
            lookup[copy1_chr] = "copy1"
        if copy2_chr:
            lookup[copy2_chr] = "copy2"
    return lookup


def add_source_copy_to_provenance(provenance_rows, corrected_homolog_rows):
    source_lookup = build_source_copy_lookup(corrected_homolog_rows)
    enriched = []
    for row in provenance_rows:
        copy_row = dict(row)
        copy_row["source_copy"] = source_lookup.get(copy_row.get("chr", ""), "")
        enriched.append(copy_row)
    return enriched


def build_reference_query_summary(rows, out_path):
    agg = defaultdict(lambda: {"aln_len": 0, "matches": 0, "covered": []})
    for row in rows:
        key = (row["target"], row["query"])
        agg[key]["aln_len"] += row["aln_len"]
        agg[key]["matches"] += row["matches"]
        agg[key]["covered"].append((row["target_start"], row["target_end"]))

    out_rows = []
    for (ref, query), stats in sorted(agg.items(), key=lambda kv: (natural_key(kv[0][0]), natural_key(kv[0][1]))):
        out_rows.append({
            "reference_chr": ref,
            "query_chr": query,
            "aligned_block_bp": str(stats["aln_len"]),
            "matching_bp": str(stats["matches"]),
            "aligned_covered_bp": str(interval_union_bp(stats["covered"])),
        })
    write_tsv(out_path, out_rows, ["reference_chr", "query_chr", "aligned_block_bp", "matching_bp", "aligned_covered_bp"])


def write_hap_folder_candidate_tables(out_dir: Path, candidate_rows):
    hap_fields = ["candidate_name", "quartet_id", "source_chr", "notes"]
    hap_rows = {"hapA1": [], "hapA2": [], "hapB1": [], "hapB2": []}
    key_map = {
        "hapA1": "candidate_A1",
        "hapA2": "candidate_A2",
        "hapB1": "candidate_B1",
        "hapB2": "candidate_B2",
    }
    for idx, row in enumerate(sorted(candidate_rows, key=lambda r: natural_key(r["quartet_id"])), start=1):
        for hap in ["hapA1", "hapA2", "hapB1", "hapB2"]:
            source_chr = row[key_map[hap]]
            hap_rows[hap].append({
                "candidate_name": f"{hap}_{idx}",
                "quartet_id": row["quartet_id"],
                "source_chr": source_chr,
                "notes": "candidate only; freeze after manual dotplot/heatmap review",
            })
    for hap, rows in hap_rows.items():
        hap_dir = out_dir / hap
        hap_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(hap_dir / "candidate_members.tsv", rows, hap_fields)
        (hap_dir / "README.txt").write_text(
            "\n".join([
                f"{hap} candidate chromosome table",
                "",
                "`candidate_members.tsv` stores chromosome-level candidate names such as",
                f"`{hap}_1`, `{hap}_2`, ..., not subfolder or FASTA names.",
                "The workflow remains candidate-only and does not rewrite FASTA.",
                "",
            ]),
            encoding="utf-8",
        )


def validate_quartets(quartets, quartet_map):
    warnings = []
    if len(quartets) != 8:
        warnings.append(f"Expected 8 homologous quartets from an n=8 reference, observed {len(quartets)}.")
    for quartet in quartets:
        members = quartet_map.get(quartet["quartet_id"], [])
        if len(members) != 4:
            warnings.append(f"{quartet['quartet_id']} has {len(members)} pzi members, expected 4.")
    return warnings


def parse_args():
    p = argparse.ArgumentParser(description="Build haplotype phase-candidate tables for A1/A2/B1/B2.")
    p.add_argument("--root-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pzi-chrom-fasta", required=True)
    p.add_argument("--pcam-fasta", required=True)
    p.add_argument("--pcam-assembly-report", required=True)
    p.add_argument("--ppa-fasta", required=True)
    p.add_argument("--ppa-assembly-report", required=True)
    p.add_argument("--chr-hap-provenance", required=True)
    p.add_argument("--corrected-homolog-pairs", required=True)
    p.add_argument("--pzi-vs-pcam-paf", required=True)
    p.add_argument("--pzi-vs-ppa-paf", required=True)
    p.add_argument("--n8-reference", action="append", default=[], help="label|fasta|assembly_report|paf or label|fasta|assembly_report|paf|expected_count")
    p.add_argument("--outgroup-reference", action="append", default=[], help="label|fasta|assembly_report|paf or label|fasta|assembly_report|paf|expected_count")
    p.add_argument("--focus-chrs", default="", help="Optional comma-separated PZI chromosomes for extra outgroup segment tables. Empty means all chr01-32.")
    p.add_argument("--segment-bin-size", type=int, default=500000)
    p.add_argument("--chr32-self-paf", required=True)
    p.add_argument("--min-mapq", type=int, default=10)
    p.add_argument("--min-aln-len", type=int, default=10000)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chr_fasta = Path(args.pzi_chrom_fasta)
    pcam_fasta = Path(args.pcam_fasta)
    ppa_fasta = Path(args.ppa_fasta)
    provenance_rows = read_tsv(Path(args.chr_hap_provenance))
    corrected_homolog_rows = read_tsv(Path(args.corrected_homolog_pairs))
    provenance_rows = add_source_copy_to_provenance(provenance_rows, corrected_homolog_rows)

    # PAF-derived summaries.
    pcam_rows = parse_paf(Path(args.pzi_vs_pcam_paf), "pcam", args.min_mapq, args.min_aln_len)
    ppa_rows = parse_paf(Path(args.pzi_vs_ppa_paf), "ppa", args.min_mapq, args.min_aln_len)
    chr32_rows = parse_paf(Path(args.chr32_self_paf), "self", args.min_mapq, args.min_aln_len)

    pcam_chr_headers, pcam_manifest = main_reference_order(
        read_fasta_headers(Path(args.pcam_fasta)),
        Path(args.pcam_assembly_report),
        8,
    )
    ppa_chr_headers, ppa_manifest = main_reference_order(
        read_fasta_headers(Path(args.ppa_fasta)),
        Path(args.ppa_assembly_report),
        16,
    )
    pzi_chr_headers = [h.split()[0] for h in read_fasta_headers(chr_fasta)]
    chr_order = [c for c in pzi_chr_headers if c.startswith("chr")]
    write_tsv(out_dir / "pcam_main_chromosomes.tsv", pcam_manifest, [
        "fasta_chr", "sequence_name", "assigned_molecule", "genbank_accn", "refseq_accn", "length_bp", "assembly_unit"
    ])
    write_tsv(out_dir / "ppa_main_chromosomes.tsv", ppa_manifest, [
        "fasta_chr", "sequence_name", "assigned_molecule", "genbank_accn", "refseq_accn", "length_bp", "assembly_unit"
    ])

    pcam_summary, pcam_best = summarize_pzi_vs_pcam(pcam_rows, pzi_chr_headers, pcam_chr_headers)
    write_tsv(out_dir / "pzi_vs_pcam_chr_alignment.tsv", pcam_summary, [
        "pzi_chr", "pcam_chr", "pzi_length_bp", "pcam_length_bp", "aligned_block_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
    ])
    write_tsv(out_dir / "pzi_vs_pcam_best_hits.tsv", pcam_best, [
        "pzi_chr", "pcam_chr", "pzi_length_bp", "pcam_length_bp", "aligned_block_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
    ])

    ppa_summary, ppa_best = summarize_pzi_vs_ppa(ppa_rows, pzi_chr_headers, ppa_chr_headers)
    write_tsv(out_dir / "pzi_vs_ppa_chr_alignment.tsv", ppa_summary, [
        "pzi_chr", "ppa_chr", "pzi_length_bp", "ppa_length_bp", "aligned_block_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
    ])
    write_tsv(out_dir / "pzi_vs_ppa_best_hits.tsv", ppa_best, [
        "pzi_chr", "ppa_chr", "pzi_length_bp", "ppa_length_bp", "aligned_block_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_pzi_chr"
    ])

    quartet_map = {}
    quartets, quartet_members = build_homolog_quartets(pcam_summary, provenance_rows, quartet_map)
    write_tsv(out_dir / "homologous_quartets.tsv", quartets, [
        "quartet_id", "pcam_chr", "pzi_chr_1", "pzi_chr_2", "pzi_chr_3", "pzi_chr_4"
    ])
    quartet_warnings = validate_quartets(quartets, quartet_map)
    quartet_check_rows = []
    for quartet in quartets:
        members = quartet_map.get(quartet["quartet_id"], [])
        quartet_check_rows.append({
            "quartet_id": quartet["quartet_id"],
            "pcam_chr": quartet["pcam_chr"],
            "member_count": str(len(members)),
            "check_status": "pass" if len(members) == 4 else "review",
            "members": ",".join(members),
        })
    write_tsv(out_dir / "homologous_quartets.check.tsv", quartet_check_rows, [
        "quartet_id", "pcam_chr", "member_count", "check_status", "members"
    ])
    (out_dir / "homologous_quartets.warnings.txt").write_text(
        "\n".join(quartet_warnings) + ("\n" if quartet_warnings else "No quartet-count warnings.\n"),
        encoding="utf-8",
    )
    write_tsv(out_dir / "quartet_members_with_copy.tsv", quartet_members, [
        "quartet_id", "quartet_rank", "pcam_chr", "pzi_chr", "source_copy", "dominant_hap", "copy_provenance_status", "copy_dominance", "copy_top_query"
    ])

    pcam_review = build_reference_review(pcam_best, "pcam")
    ppa_review = build_reference_review(ppa_best, "ppa")
    write_tsv(out_dir / "pcam_rearrangement_review.tsv", pcam_review, [
        "reference_species", "reference_chr", "pzi_members", "member_count", "reference_length_bp",
        "pzi_member_total_bp", "length_ratio_to_reference", "event_class", "eej_ncf_rct_candidate",
        "review_status", "review_note"
    ])
    write_tsv(out_dir / "ppa_rearrangement_review.tsv", ppa_review, [
        "reference_species", "reference_chr", "pzi_members", "member_count", "reference_length_bp",
        "pzi_member_total_bp", "length_ratio_to_reference", "event_class", "eej_ncf_rct_candidate",
        "review_status", "review_note"
    ])
    write_tsv(out_dir / "trispecies_pzi_reference_best_hits.tsv", build_trispecies_pzi_table(chr_order, pcam_best, ppa_best), [
        "pzi_chr", "pcam_chr", "pcam_coverage_ratio", "pcam_orientation", "ppa_chr",
        "ppa_coverage_ratio", "ppa_orientation", "review_status", "review_note"
    ])

    n8_specs = [parse_reference_spec(text) for text in args.n8_reference]
    outgroup_specs = [parse_reference_spec(text) for text in args.outgroup_reference]
    focus_chrs = parse_chr_list(args.focus_chrs) or chr_order
    n8_outputs = build_n8_reference_outputs(
        args,
        out_dir,
        n8_specs,
        pzi_chr_headers,
        chr_order,
        args.min_mapq,
        args.min_aln_len,
    ) if n8_specs else {}
    outgroup_outputs = build_outgroup_outputs(
        out_dir,
        outgroup_specs,
        pzi_chr_headers,
        focus_chrs,
        args.min_mapq,
        args.min_aln_len,
        args.segment_bin_size,
    ) if outgroup_specs else {}
    consensus_rows = build_multi_reference_consensus(chr_order, ppa_best, pcam_best, n8_outputs)
    write_tsv(out_dir / "multi_reference_consensus_by_ppa_group.tsv", consensus_rows, [
        "ppa_chr", "pzi_members", "member_count", "pcam_vote", "n8_consensus_vote",
        "n8_consensus_count", "all_n8_votes", "status", "review_note"
    ])
    consensus_quartets, consensus_members = build_consensus_quartets(consensus_rows)
    write_tsv(out_dir / "multi_reference_consensus_quartets.tsv", consensus_quartets, [
        "consensus_quartet_id", "consensus_key", "ppa_groups", "member_count",
        "pzi_chr_1", "pzi_chr_2", "pzi_chr_3", "pzi_chr_4", "status", "review_note"
    ])
    write_tsv(out_dir / "multi_reference_consensus_members.tsv", consensus_members, [
        "consensus_quartet_id", "consensus_key", "rank", "pzi_chr", "ppa_groups", "status"
    ])

    consensus_color_by_pzi = color_by_pzi_from_consensus(consensus_quartets)
    color_by_pzi, color_by_pcam = build_color_maps_from_quartets(quartets, quartet_map)
    if consensus_color_by_pzi:
        color_by_pzi.update(consensus_color_by_pzi)
    color_by_ppa = build_reference_color_map(ppa_best, "ppa", color_by_pzi)
    write_axis_colors(out_dir / "pzi_vs_pcam_axis_colors.tsv", color_by_pzi, color_by_pcam, "pcam")
    write_axis_colors(out_dir / "pzi_vs_ppa_axis_colors.tsv", color_by_pzi, color_by_ppa, "ppa")
    write_axis_segments(out_dir / "pzi_vs_pcam_axis_segments.tsv", pcam_rows, color_by_pzi, color_by_pcam)
    write_axis_segments(out_dir / "pzi_vs_ppa_axis_segments.tsv", ppa_rows, color_by_pzi, color_by_ppa)
    for label, data in n8_outputs.items():
        color_by_ref = build_reference_color_map(data["best"], label, color_by_pzi)
        write_axis_colors(out_dir / f"pzi_vs_{label}_axis_colors.tsv", color_by_pzi, color_by_ref, label)
        write_axis_segments(out_dir / f"pzi_vs_{label}_axis_segments.tsv", data["paf_rows"], color_by_pzi, color_by_ref)
    for label, data in outgroup_outputs.items():
        color_by_ref = build_reference_color_map(data["best"], label, color_by_pzi)
        write_axis_colors(out_dir / f"pzi_vs_{label}_axis_colors.tsv", color_by_pzi, color_by_ref, label)
        write_axis_segments(out_dir / f"pzi_vs_{label}_axis_segments.tsv", data["paf_rows"], color_by_pzi, color_by_ref)

    provenance_lookup = {row["chr"]: row for row in provenance_rows}
    candidate_rows = infer_haplotype_candidates(quartet_map, provenance_rows, corrected_homolog_rows)
    write_tsv(out_dir / "candidate_haplotype_naming.tsv", candidate_rows, [
        "quartet_id", "member_1", "member_2", "member_3", "member_4",
        "member_1_hap", "member_2_hap", "member_3_hap", "member_4_hap",
        "member_1_copy", "member_2_copy", "member_3_copy", "member_4_copy",
        "candidate_A1", "candidate_A2", "candidate_B1", "candidate_B2",
        "candidate_A1_copy", "candidate_A2_copy", "candidate_B1_copy", "candidate_B2_copy",
        "candidate_basis", "subgenome_label_warning", "status",
    ])
    write_hap_folder_candidate_tables(out_dir, candidate_rows)

    # copy1/copy2 overlay at quartet level.
    overlay_rows = []
    for row in quartet_members:
        chr_name = row["pzi_chr"]
        prov = provenance_lookup.get(chr_name, {})
        overlay_rows.append({
            **row,
            "provenance_dominant_hap": prov.get("dominant_hap", ""),
            "dominance": prov.get("dominance", ""),
            "provenance_status": prov.get("status", ""),
        })
    write_tsv(out_dir / "quartet_overlay_with_copy_provenance.tsv", overlay_rows, [
        "quartet_id", "quartet_rank", "pcam_chr", "pzi_chr", "source_copy", "dominant_hap", "copy_provenance_status", "copy_dominance", "copy_top_query",
        "provenance_dominant_hap", "dominance", "provenance_status"
    ])

    # chr32 self-synteny matrix.
    self_matrix = build_all_vs_all_matrix(chr32_rows, chr_order)
    write_tsv(out_dir / "chr32_all_vs_all_matrix.tsv", self_matrix, ["query_chr", "target_chr", "aligned_bp"])

    # Heatmap-friendly summaries.
    build_reference_query_summary(pcam_rows, out_dir / "pzi_vs_pcam_alignment_summary.tsv")
    build_reference_query_summary(ppa_rows, out_dir / "pzi_vs_ppa_alignment_summary.tsv")
    build_reference_query_summary(chr32_rows, out_dir / "chr32_self_alignment_summary.tsv")

    # Quartet statistics and JSON manifest for downstream manual review.
    quartet_meta = []
    for quartet_id, members in quartet_map.items():
        quartet_meta.append({
            "quartet_id": quartet_id,
            "member_count": str(len(members)),
            "members": ",".join(members),
        })
    write_tsv(out_dir / "quartet_manifest.tsv", quartet_meta, ["quartet_id", "member_count", "members"])
    (out_dir / "run_inputs.json").write_text(json.dumps({
        "root_dir": args.root_dir,
        "pzi_chrom_fasta": args.pzi_chrom_fasta,
        "pcam_fasta": args.pcam_fasta,
        "pcam_assembly_report": args.pcam_assembly_report,
        "ppa_fasta": args.ppa_fasta,
        "ppa_assembly_report": args.ppa_assembly_report,
        "chr_hap_provenance": args.chr_hap_provenance,
        "pzi_vs_pcam_paf": args.pzi_vs_pcam_paf,
        "pzi_vs_ppa_paf": args.pzi_vs_ppa_paf,
        "n8_references": [spec["label"] for spec in n8_specs],
        "outgroup_references": [spec["label"] for spec in outgroup_specs],
        "focus_chrs": focus_chrs,
        "segment_bin_size": args.segment_bin_size,
        "chr32_self_paf": args.chr32_self_paf,
        "pcam_main_chromosomes": len(pcam_chr_headers),
        "ppa_main_chromosomes": len(ppa_chr_headers),
        "provenance_rows": len(provenance_rows),
        "corrected_homolog_rows": len(corrected_homolog_rows),
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Wrote haploid phase-candidate tables into {out_dir}")


if __name__ == "__main__":
    main()
