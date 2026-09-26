#!/usr/bin/env python3
"""Catalogue, compare and validate large PZI tandem-repeat arrays.

This workflow makes deliberately limited claims.  Flank-pair mappings estimate
the distance between two putatively orthologous flanks; they do not by
themselves prove the evolutionary origin of a repeat.  ONT evidence tests
assembly support at the two array junctions, while a complete spanning read is
reported only when one actually covers the entire selected array.
"""

import argparse
import bisect
import collections
import csv
import gzip
import json
import subprocess
from pathlib import Path


PZI_TARGETS = ("pziA1", "pziA2", "pziB1", "pziB2")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(a, b) for a, b in result]


def span(intervals):
    return sum(end - start for start, end in intervals)


def intersection(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def fai_records(path):
    result = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                result[fields[0]] = tuple(map(int, fields[1:5]))
    if not result:
        raise ValueError(f"empty or invalid FASTA index: {path}")
    return result


def fetch_fasta(fasta, index, seqid, start, end):
    """Fetch a zero-based half-open interval from an uncompressed indexed FASTA."""
    if seqid not in index:
        raise KeyError(f"{seqid} absent from {fasta}.fai")
    length, byte_offset, line_bases, line_width = index[seqid]
    start, end = max(0, start), min(length, end)
    if end <= start:
        return ""
    pieces = []
    with open(fasta, "rb") as handle:
        position = start
        while position < end:
            within_line = position % line_bases
            take = min(end - position, line_bases - within_line)
            offset = byte_offset + (position // line_bases) * line_width + within_line
            handle.seek(offset)
            pieces.append(handle.read(take).decode("ascii"))
            position += take
    return "".join(pieces).upper()


def trf_calls(path, min_score):
    by_seq = collections.defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                start, end = int(row["start"]), int(row["end"])
                score, period = int(row["score"]), int(row["period_bp"])
            except (KeyError, ValueError):
                continue
            if end > start and score >= min_score:
                by_seq[row["seqid"]].append({
                    "start": start, "end": end, "score": score, "period_bp": period,
                    "motif": row.get("consensus_motif", ""),
                })
    return by_seq


def make_arrays(calls, join_gap):
    arrays = []
    for seqid, values in calls.items():
        current = []
        end = None
        for row in sorted(values, key=lambda item: (item["start"], item["end"])):
            if end is None or row["start"] <= end + join_gap:
                current.append(row)
                end = max(end if end is not None else row["end"], row["end"])
                continue
            arrays.append((seqid, current))
            current, end = [row], row["end"]
        if current:
            arrays.append((seqid, current))
    result = []
    for seqid, items in arrays:
        intervals = merge((row["start"], row["end"]) for row in items)
        longest = max(items, key=lambda row: (row["end"] - row["start"], row["score"]))
        result.append({
            "seqid": seqid,
            "start": min(row["start"] for row in items),
            "end": max(row["end"] for row in items),
            "array_span_bp": max(row["end"] for row in items) - min(row["start"] for row in items),
            "trf_union_bp": span(intervals),
            "trf_interval_count": len(items),
            "max_trf_score": max(row["score"] for row in items),
            "representative_period_bp": longest["period_bp"],
            "representative_consensus_motif": longest["motif"],
        })
    return result


def command_catalogue(args):
    manifest = {row["target"]: row for row in read_tsv(args.mask_manifest)}
    missing = [target for target in PZI_TARGETS if target not in manifest]
    if missing:
        raise SystemExit("ERROR: missing PZI targets in mask manifest: " + ", ".join(missing))
    rows, fasta_records = [], {}
    for target in PZI_TARGETS:
        fasta = Path(manifest[target]["mask_input_fasta"])
        trf = args.tandem_root / target / f"{target}.trf.tsv"
        fai = Path(str(fasta) + ".fai")
        for path in (fasta, fai, trf):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"ERROR: missing input for {target}: {path}")
        index = fai_records(fai)
        fasta_records[target] = (fasta, index)
        arrays = make_arrays(trf_calls(trf, args.min_score), args.join_gap_bp)
        arrays = [row for row in arrays if row["trf_union_bp"] >= args.min_array_union_bp]
        arrays.sort(key=lambda row: (-row["trf_union_bp"], -row["array_span_bp"], row["seqid"], row["start"]))
        for rank, row in enumerate(arrays, 1):
            row.update({"target": target, "array_id": f"{target}_array{rank:04d}", "rank_within_target": rank})
            length = index.get(row["seqid"], (0,))[0]
            row["flank_bp"] = args.flank_bp
            row["flank_eligible"] = int(row["start"] >= args.flank_bp and row["end"] + args.flank_bp <= length)
            rows.append(row)
    fields = [
        "target", "array_id", "rank_within_target", "seqid", "start", "end", "array_span_bp",
        "trf_union_bp", "trf_interval_count", "max_trf_score", "representative_period_bp",
        "representative_consensus_motif", "flank_bp", "flank_eligible",
    ]
    write_tsv(args.output_dir / "pzi_tandem_array_catalogue.tsv", fields, rows)
    selected = [row for row in rows if row["rank_within_target"] <= args.top_arrays_per_target and row["flank_eligible"]]
    write_tsv(args.output_dir / "pzi_tandem_array_candidates.tsv", fields, selected)
    fasta_path = args.output_dir / "pzi_tandem_array_flanks.fa"
    with open(fasta_path, "w", encoding="utf-8") as handle:
        for row in selected:
            fasta, index = fasta_records[row["target"]]
            for side, start, end in (
                ("L", int(row["start"]) - args.flank_bp, int(row["start"])),
                ("R", int(row["end"]), int(row["end"]) + args.flank_bp),
            ):
                sequence = fetch_fasta(fasta, index, row["seqid"], start, end)
                if len(sequence) != args.flank_bp:
                    raise SystemExit(f"ERROR: failed to extract complete {side} flank for {row['array_id']}")
                handle.write(f">{row['target']}|{row['array_id']}|{side}|{row['seqid']}:{start}-{end}\n")
                for offset in range(0, len(sequence), 80):
                    handle.write(sequence[offset:offset + 80] + "\n")
    metadata = {
        "trf_min_score": args.min_score, "array_join_gap_bp": args.join_gap_bp,
        "min_array_union_bp": args.min_array_union_bp, "flank_bp": args.flank_bp,
        "top_arrays_per_target": args.top_arrays_per_target,
        "interpretation": "Candidates are large PZI TRF-array clusters selected for locus-level comparison, not all independently evolved repeat families.",
    }
    with open(args.output_dir / "pzi_tandem_array_catalogue_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    print(json.dumps({"catalogued_arrays": len(rows), "flank_candidates": len(selected), "output_dir": str(args.output_dir)}))


def parse_paf(path, min_mapq, min_query_coverage, min_identity, ambiguity_relative_match):
    hits = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                qlen, qs, qe = int(fields[1]), int(fields[2]), int(fields[3])
                ts, te = int(fields[7]), int(fields[8])
                matches, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
            except ValueError:
                continue
            qcov = (qe - qs) / qlen if qlen else 0
            identity = matches / block if block else 0
            if mapq < min_mapq or qcov < min_query_coverage or identity < min_identity:
                continue
            hits[fields[0]].append({
                "target_seqid": fields[5], "target_start": ts, "target_end": te,
                "strand": fields[4], "matches": matches, "block_bp": block,
                "mapq": mapq, "query_coverage": qcov, "identity": identity,
            })
    best, ambiguous = {}, set()
    for query, values in hits.items():
        values.sort(key=lambda row: (-row["mapq"], -row["matches"], -row["query_coverage"], -row["identity"]))
        best[query] = values[0]
        top = values[0]
        # A repeat-array flank must resolve to one genomic context.  A second
        # comparably strong location is intentionally treated as unresolved.
        for candidate in values[1:]:
            same_locus = (
                candidate["target_seqid"] == top["target_seqid"]
                and candidate["target_start"] == top["target_start"]
                and candidate["target_end"] == top["target_end"]
                and candidate["strand"] == top["strand"]
            )
            if not same_locus and candidate["matches"] >= top["matches"] * ambiguity_relative_match:
                ambiguous.add(query)
                break
    return best, ambiguous


def trf_unions_for_target(tandem_root, target, min_score, cache):
    if target not in cache:
        path = tandem_root / target / f"{target}.trf.tsv"
        if not path.is_file():
            cache[target] = {}
        else:
            calls = trf_calls(path, min_score)
            cache[target] = {seqid: merge((row["start"], row["end"]) for row in rows) for seqid, rows in calls.items()}
    return cache[target]


def command_summarize_flanks(args):
    candidates = read_tsv(args.candidates)
    candidate_by_key = {(row["target"], row["array_id"]): row for row in candidates}
    compare_rows = read_tsv(args.comparison_manifest)
    required = {"query_target", "target", "paf"}
    if not compare_rows or not required.issubset(compare_rows[0]):
        raise SystemExit("ERROR: comparison manifest requires query_target, target and paf columns")
    trf_cache, rows = {}, []
    for compare in compare_rows:
        query_target, target, paf = compare["query_target"], compare["target"], Path(compare["paf"])
        if not paf.is_file() or paf.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing PAF: {paf}")
        hit, ambiguous = parse_paf(
            paf, args.min_mapq, args.min_query_coverage, args.min_identity,
            args.ambiguity_relative_match,
        )
        target_trf = trf_unions_for_target(args.tandem_root, target, args.min_trf_score, trf_cache)
        for key, candidate in candidate_by_key.items():
            if key[0] != query_target:
                continue
            left_name = f"{query_target}|{candidate['array_id']}|L|{candidate['seqid']}:{int(candidate['start']) - int(candidate['flank_bp'])}-{candidate['start']}"
            right_name = f"{query_target}|{candidate['array_id']}|R|{candidate['seqid']}:{candidate['end']}-{int(candidate['end']) + int(candidate['flank_bp'])}"
            left, right = hit.get(left_name), hit.get(right_name)
            base = {
                "query_target": query_target, "array_id": candidate["array_id"], "query_seqid": candidate["seqid"],
                "query_start": candidate["start"], "query_end": candidate["end"], "query_array_span_bp": candidate["array_span_bp"],
                "query_trf_union_bp": candidate["trf_union_bp"], "target": target,
            }
            if not left or not right or left_name in ambiguous or right_name in ambiguous:
                state = "unresolved_nonunique_flank" if left_name in ambiguous or right_name in ambiguous else "unresolved_missing_qualified_flank"
                rows.append({**base, "callable": 0, "state": state, "target_seqid": "", "target_gap_bp": "", "target_gap_to_query_span": "", "target_trf_union_bp": "", "target_trf_fraction_of_gap": "", "left_mapq": left["mapq"] if left else "", "right_mapq": right["mapq"] if right else "", "left_identity": f"{left['identity']:.6f}" if left else "", "right_identity": f"{right['identity']:.6f}" if right else ""})
                continue
            if left["target_seqid"] != right["target_seqid"] or left["strand"] != right["strand"]:
                rows.append({**base, "callable": 0, "state": "unresolved_discordant_flanks", "target_seqid": "", "target_gap_bp": "", "target_gap_to_query_span": "", "target_trf_union_bp": "", "target_trf_fraction_of_gap": "", "left_mapq": left["mapq"], "right_mapq": right["mapq"], "left_identity": f"{left['identity']:.6f}", "right_identity": f"{right['identity']:.6f}"})
                continue
            if left["strand"] == "+":
                gap_start, gap_end = left["target_end"], right["target_start"]
            else:
                gap_start, gap_end = right["target_end"], left["target_start"]
            if gap_end < gap_start:
                rows.append({**base, "callable": 0, "state": "unresolved_noncollinear_flanks", "target_seqid": left["target_seqid"], "target_gap_bp": "", "target_gap_to_query_span": "", "target_trf_union_bp": "", "target_trf_fraction_of_gap": "", "left_mapq": left["mapq"], "right_mapq": right["mapq"], "left_identity": f"{left['identity']:.6f}", "right_identity": f"{right['identity']:.6f}"})
                continue
            target_gap = gap_end - gap_start
            target_trf_bp = intersection([(gap_start, gap_end)], target_trf.get(left["target_seqid"], []))
            ratio = target_gap / int(candidate["array_span_bp"]) if int(candidate["array_span_bp"]) else 0
            state = "contracted_target_interval" if ratio <= args.contracted_ratio else ("expanded_target_interval" if ratio >= args.expanded_ratio else "similar_sized_target_interval")
            rows.append({**base, "callable": 1, "state": state, "target_seqid": left["target_seqid"], "target_gap_bp": target_gap, "target_gap_to_query_span": f"{ratio:.8f}", "target_trf_union_bp": target_trf_bp, "target_trf_fraction_of_gap": f"{target_trf_bp / target_gap:.8f}" if target_gap else "0", "left_mapq": left["mapq"], "right_mapq": right["mapq"], "left_identity": f"{left['identity']:.6f}", "right_identity": f"{right['identity']:.6f}"})

    fields = [
        "query_target", "array_id", "query_seqid", "query_start", "query_end", "query_array_span_bp", "query_trf_union_bp",
        "target", "callable", "state", "target_seqid", "target_gap_bp", "target_gap_to_query_span", "target_trf_union_bp",
        "target_trf_fraction_of_gap", "left_mapq", "right_mapq", "left_identity", "right_identity",
    ]
    write_tsv(args.output_dir / "pzi_tandem_flank_comparisons.tsv", fields, rows)
    by_array = collections.defaultdict(list)
    for row in rows:
        by_array[(row["query_target"], row["array_id"])].append(row)
    summary = []
    for (query, array_id), values in sorted(by_array.items()):
        pzi = [r for r in values if r["target"].startswith("pzi")]
        psp = [r for r in values if r["target"].startswith("psp")]
        outgroups = [r for r in values if r["target"] in ("ppaA", "ppaB", "pcam", "ppe")]
        pzi_callable = [r for r in pzi if int(r["callable"])]
        psp_callable = [r for r in psp if int(r["callable"])]
        outgroup_callable = [r for r in outgroups if int(r["callable"])]
        pzi_similar = [r for r in pzi_callable if r["state"] == "similar_sized_target_interval"]
        psp_contracted = [r for r in psp_callable if r["state"] == "contracted_target_interval"]
        outgroup_contracted = [r for r in outgroup_callable if r["state"] == "contracted_target_interval"]
        if pzi_similar and psp_contracted and outgroup_contracted:
            conclusion = "PZI_homologous_array_with_PSP_and_outgroup_contracted_interval"
        elif pzi_similar and psp_contracted:
            conclusion = "PZI_homologous_array_with_PSP_contracted_interval"
        elif pzi_callable or psp_callable:
            conclusion = "callable_but_not_diagnostic"
        else:
            conclusion = "unresolved"
        summary.append({
            "query_target": query, "array_id": array_id,
            "pzi_counterpart_callable": len(pzi_callable), "pzi_counterpart_similar_sized": len(pzi_similar),
            "psp_counterpart_callable": len(psp_callable), "psp_counterpart_contracted": len(psp_contracted),
            "outgroup_counterpart_callable": len(outgroup_callable), "outgroup_counterpart_contracted": len(outgroup_contracted),
            "locus_conclusion": conclusion,
        })
    write_tsv(args.output_dir / "pzi_tandem_locus_summary.tsv", list(summary[0]) if summary else [], summary)
    print(json.dumps({"comparison_rows": len(rows), "array_summaries": len(summary), "output_dir": str(args.output_dir)}))


def cigar_span(cigar):
    digits, reference = "", 0
    for char in cigar:
        if char.isdigit():
            digits += char
            continue
        if not digits:
            return 0
        value = int(digits)
        if char in "MDN=X":
            reference += value
        digits = ""
    return reference if not digits else 0


def command_ont(args):
    candidates = [row for row in read_tsv(args.candidates) if row["target"] == args.target and int(row.get("flank_eligible", 0))]
    candidates = sorted(candidates, key=lambda row: int(row["rank_within_target"]))[:args.top_arrays]
    if not candidates:
        raise SystemExit(f"ERROR: no flank-eligible candidates for {args.target}")
    bam = Path(args.bam)
    if not bam.is_file() or bam.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing BAM: {bam}")
    boundaries = collections.defaultdict(list)
    evidence = {}
    for index, row in enumerate(candidates):
        seqid, start, end = row["seqid"], int(row["start"]), int(row["end"])
        boundaries[seqid].append((start, index, "left"))
        boundaries[seqid].append((end, index, "right"))
        evidence[index] = {"left": set(), "right": set(), "whole": set()}
    positions = {seqid: [item[0] for item in sorted(values)] for seqid, values in boundaries.items()}
    boundaries = {seqid: sorted(values) for seqid, values in boundaries.items()}
    command = ["samtools", "view", "-F", str(args.exclude_flags), "-q", str(args.min_mapq), str(bam)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6 or fields[2] not in boundaries:
            continue
        start = int(fields[3]) - 1
        end = start + cigar_span(fields[5])
        if end <= start:
            continue
        seqid = fields[2]
        lo = bisect.bisect_left(positions[seqid], start + args.flank_bp)
        hi = bisect.bisect_right(positions[seqid], end - args.flank_bp)
        for boundary, index, side in boundaries[seqid][lo:hi]:
            evidence[index][side].add(fields[0])
        for index, row in enumerate(candidates):
            if row["seqid"] != seqid:
                continue
            arr_start, arr_end = int(row["start"]), int(row["end"])
            if start <= arr_start - args.flank_bp and end >= arr_end + args.flank_bp:
                evidence[index]["whole"].add(fields[0])
    if process.wait() != 0:
        raise SystemExit(f"ERROR: samtools view failed for {bam}")
    rows = []
    for index, row in enumerate(candidates):
        left, right, whole = evidence[index]["left"], evidence[index]["right"], evidence[index]["whole"]
        both = left & right
        junction_pass = len(left) >= args.min_junction_reads and len(right) >= args.min_junction_reads
        rows.append({
            "target": args.target, "array_id": row["array_id"], "seqid": row["seqid"], "start": row["start"], "end": row["end"],
            "array_span_bp": row["array_span_bp"], "trf_union_bp": row["trf_union_bp"], "rank_within_target": row["rank_within_target"],
            "left_junction_spanning_primary_reads": len(left), "right_junction_spanning_primary_reads": len(right),
            "same_read_supporting_both_junctions": len(both), "complete_array_spanning_primary_reads": len(whole),
            "two_junctions_supported": int(junction_pass),
            "support_interpretation": "two_array_junctions_supported" if junction_pass else "insufficient_unique_primary_junction_support",
            "caveat": "Shared tetraploid ONT reads and repetitive array interiors can cross-map or have low MAPQ; junction support validates assembly adjacency, not exact repeat copy number.",
        })
    fields = list(rows[0])
    write_tsv(args.output, fields, rows)
    print(json.dumps({"target": args.target, "arrays_tested": len(rows), "two_junctions_supported": sum(int(r["two_junctions_supported"]) for r in rows), "output": str(args.output)}))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    catalogue = sub.add_parser("catalogue", help="Cluster PZI TRF calls into large arrays and extract flanks.")
    catalogue.add_argument("--mask-manifest", type=Path, required=True)
    catalogue.add_argument("--tandem-root", type=Path, required=True)
    catalogue.add_argument("--output-dir", type=Path, required=True)
    catalogue.add_argument("--min-score", type=int, default=50)
    catalogue.add_argument("--join-gap-bp", type=int, default=10000)
    catalogue.add_argument("--min-array-union-bp", type=int, default=100000)
    catalogue.add_argument("--top-arrays-per-target", type=int, default=25)
    catalogue.add_argument("--flank-bp", type=int, default=5000)
    catalogue.set_defaults(func=command_catalogue)
    summary = sub.add_parser("summarize-flanks", help="Classify flank-pair mappings at candidate array loci.")
    summary.add_argument("--candidates", type=Path, required=True)
    summary.add_argument("--comparison-manifest", type=Path, required=True)
    summary.add_argument("--tandem-root", type=Path, required=True)
    summary.add_argument("--output-dir", type=Path, required=True)
    summary.add_argument("--min-trf-score", type=int, default=50)
    summary.add_argument("--min-mapq", type=int, default=20)
    summary.add_argument("--min-query-coverage", type=float, default=0.80)
    summary.add_argument("--min-identity", type=float, default=0.70)
    summary.add_argument("--ambiguity-relative-match", type=float, default=0.95,
                         help="A second qualified flank hit at this fraction of the best match is unresolved.")
    summary.add_argument("--contracted-ratio", type=float, default=0.50)
    summary.add_argument("--expanded-ratio", type=float, default=1.50)
    summary.set_defaults(func=command_summarize_flanks)
    ont = sub.add_parser("ont", help="Count primary ONT reads supporting array junctions.")
    ont.add_argument("--candidates", type=Path, required=True)
    ont.add_argument("--target", choices=PZI_TARGETS, required=True)
    ont.add_argument("--bam", type=Path, required=True)
    ont.add_argument("--output", type=Path, required=True)
    ont.add_argument("--top-arrays", type=int, default=15)
    ont.add_argument("--flank-bp", type=int, default=1000)
    ont.add_argument("--min-mapq", type=int, default=20)
    ont.add_argument("--min-junction-reads", type=int, default=3)
    ont.add_argument("--exclude-flags", type=int, default=2308)
    ont.set_defaults(func=command_ont)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
