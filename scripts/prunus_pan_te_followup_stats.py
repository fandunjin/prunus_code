#!/usr/bin/env python3
"""Prepare TE, tandem-repeat and orthosite follow-up tables.

This script only performs descriptive, reproducible summaries.  It does not
treat phased haplotypes as independent biological replicates and it keeps the
TE/TRF overlap as an annotation relationship rather than an additive signal.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


TARGETS = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
GROUP = {t: ("PZI" if t.startswith("pzi") else "PSP" if t.startswith("psp")
             else "PPA" if t.startswith("ppa") else t.upper()) for t in TARGETS}


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def infer_mu(path: Path):
    """Read posterior mean and HPD interval from an MCMCtree .out file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    pattern = re.compile(
        r"^\s*mu\s+([0-9.eE+-]+)\s+\(\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)\)\s+\(\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)\)",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if not match:
        raise ValueError(f"could not find posterior mu in {path}")
    values = [float(x) for x in match.groups()]
    return {"source": path.name, "mu_mean": values[0], "mu_equal_tail_low": values[1],
            "mu_equal_tail_high": values[2], "mu_hpd_low": values[3], "mu_hpd_high": values[4]}


def age_from_divergence(divergence_pct, mu, mu_time_unit_mya=100.0):
    """Approximate age from pairwise divergence and a MCMCtree rate.

    MCMCtree rates here are interpreted per ``mu_time_unit_mya`` (the tree
    calibration uses 100-Myr-like time units).  RepeatMasker divergence from a
    consensus is only an approximate age proxy; intact-LTR values are the
    preferred estimate.
    """
    if mu <= 0:
        return ""
    return (float(divergence_pct) / 100.0) / (2.0 * mu / mu_time_unit_mya)


def summarize_family(args, out):
    rows = read_tsv(args.drivers)
    keep = [r for r in rows if fnum(r.get("positive_excess_bp")) > 0]
    keep.sort(key=lambda r: fnum(r.get("positive_excess_bp")), reverse=True)
    top = []
    for rank, row in enumerate(keep[: args.top_families], 1):
        top.append({
            "rank": rank,
            "class_family": row.get("class_family", ""),
            "pan_family": row.get("pan_family", ""),
            "pzi_mean_bp": row.get("pzi_mean_bp", "0"),
            "psp_mean_bp": row.get("psp_mean_bp", "0"),
            "pzi_psp_excess_bp": row.get("positive_excess_bp", "0"),
            "ppa_pcam_ppe_median_bp": row.get("ppa_pcam_ppe_median_bp", "0"),
            "expanded_in_both_subgenomes": row.get("expanded_in_both_subgenomes", "0"),
            "consistent_across_pzi_haplotypes": row.get("consistent_across_pzi_haplotypes", "0"),
            "cumulative_positive_excess_fraction": row.get("cumulative_positive_excess_fraction", ""),
            "library_low_complexity_flag": row.get("library_low_complexity_flag", ""),
        })
    write_tsv(out / "te_family_driver_top.tsv", top)

    # A compact class-level contribution table for the main figure.
    class_rows = defaultdict(lambda: {"family_count": 0, "excess_bp": 0.0})
    for row in keep:
        key = row.get("class_family", "Other")
        class_rows[key]["family_count"] += 1
        class_rows[key]["excess_bp"] += fnum(row.get("positive_excess_bp"))
    class_out = [
        {"class_family": key, "family_count": value["family_count"], "positive_excess_bp": round(value["excess_bp"])}
        for key, value in class_rows.items()
    ]
    class_out.sort(key=lambda r: r["positive_excess_bp"], reverse=True)
    write_tsv(out / "te_driver_class_summary.tsv", class_out)


def summarize_age(args, out):
    calibrations = []
    if args.species_mcmctree:
        calibrations.append({"calibration": "species14", **infer_mu(args.species_mcmctree)})
    if args.subgenome_mcmctree:
        calibrations.append({"calibration": "subgenome17", **infer_mu(args.subgenome_mcmctree)})
    if not calibrations:
        calibrations = [{"calibration": "user", "source": "--mu", "mu_mean": args.mu,
                         "mu_equal_tail_low": "", "mu_equal_tail_high": "", "mu_hpd_low": "", "mu_hpd_high": ""}]
    write_tsv(out / "clock_calibration.tsv", calibrations)
    by_cal = {row["calibration"]: row for row in calibrations}
    species = by_cal.get("species14", calibrations[0])
    subgenome = by_cal.get("subgenome17", species)
    mu = args.mu if args.mu is not None else fnum(species["mu_mean"])
    if mu <= 0:
        raise ValueError("a positive --mu or MCMCtree output is required")

    def age_fields(div):
        result = {
            "estimated_age_mya": f"{age_from_divergence(div, mu, args.mu_time_unit_mya):.6f}",
            "species_mu_age_low_mya": f"{age_from_divergence(div, fnum(species.get('mu_hpd_high'), mu), args.mu_time_unit_mya):.6f}",
            "species_mu_age_high_mya": f"{age_from_divergence(div, fnum(species.get('mu_hpd_low'), mu), args.mu_time_unit_mya):.6f}",
            "subgenome_mu_age_mya": f"{age_from_divergence(div, fnum(subgenome.get('mu_mean'), mu), args.mu_time_unit_mya):.6f}",
            "subgenome_mu_age_low_mya": f"{age_from_divergence(div, fnum(subgenome.get('mu_hpd_high'), mu), args.mu_time_unit_mya):.6f}",
            "subgenome_mu_age_high_mya": f"{age_from_divergence(div, fnum(subgenome.get('mu_hpd_low'), mu), args.mu_time_unit_mya):.6f}",
        }
        return result

    divergence = read_tsv(args.divergence)
    div_out = []
    for row in divergence:
        div = fnum(row.get("divergence_bin_mid_pct"))
        if div > args.max_divergence:
            continue
        div_out.append({
            "target": row.get("target", ""), "group": row.get("group", ""),
            "class_family": row.get("class_family", ""), "pan_family": row.get("pan_family", ""),
            "divergence_pct": f"{div:.4f}", "exclusive_bp": row.get("exclusive_bp", "0"),
            **age_fields(div),
            "age_model": "RepeatMasker divergence from consensus; approximate",
        })
    write_tsv(out / "te_divergence_age.tsv", div_out)

    intact = read_tsv(args.intact_ltr)
    intact_out = []
    for row in intact:
        div = fnum(row.get("ltr_divergence")) * 100.0
        intact_out.append({
            "target": row.get("target", ""), "group": row.get("group", ""),
            "seqid": row.get("seqid", ""), "start": row.get("start", ""), "end": row.get("end", ""),
            "length_bp": row.get("length_bp", "0"), "pan_family": row.get("pan_family", ""),
            "class_family": row.get("pan_class_family", row.get("original_class_family", "")),
            "ltr_divergence_pct": f"{div:.4f}",
            **age_fields(div),
            "ltr_identity": row.get("ltr_identity", ""),
        })
    write_tsv(out / "intact_ltr_age.tsv", intact_out)


def canonical_motif(value):
    value = (value or "").split(":", 1)[-1].upper()
    value = re.sub(r"[^ACGT]", "", value)
    if not value:
        return ""
    comp = str.maketrans("ACGT", "TGCA")
    rc = value.translate(comp)[::-1]
    return min(value, rc)


def motif_identity(a, b):
    if not a or not b:
        return 0.0
    if len(a) != len(b):
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


def summarize_motifs(args, out):
    rows = read_tsv(args.arrays)
    motifs = defaultdict(lambda: {"array_count": 0, "bp": 0, "targets": set(), "period_bp": 0})
    for row in rows:
        motif = canonical_motif(row.get("canonical_motif") or row.get("raw_consensus_motif"))
        if not motif:
            continue
        key = (int(fnum(row.get("period_bp"))), motif)
        record = motifs[key]
        record["array_count"] += 1
        record["bp"] += int(fnum(row.get("array_length_bp")))
        record["targets"].add(row.get("target", ""))
        record["period_bp"] = int(fnum(row.get("period_bp")))
    representatives = sorted(motifs.items(), key=lambda kv: kv[1]["bp"], reverse=True)[: args.max_motifs]
    # Greedy clustering among equal-period motifs.  This is deliberately
    # conservative and is a screening cluster, not a phylogenetic inference.
    clusters = []
    for (period, motif), record in representatives:
        assigned = None
        for cluster in clusters:
            if abs(period - cluster["period_bp"]) > max(2, int(0.05 * max(period, cluster["period_bp"]))):
                continue
            if motif_identity(motif, cluster["representative_motif"]) >= args.motif_identity:
                assigned = cluster
                break
        if assigned is None:
            assigned = {"cluster_id": f"TRFCL_{len(clusters) + 1:04d}", "period_bp": period,
                        "representative_motif": motif, "bp": 0, "array_count": 0,
                        "targets": set(), "member_count": 0, "members": []}
            clusters.append(assigned)
        assigned["bp"] += record["bp"]
        assigned["array_count"] += record["array_count"]
        assigned["targets"].update(record["targets"])
        assigned["member_count"] += 1
        assigned["members"].append({"motif": motif, "bp": record["bp"], "array_count": record["array_count"]})
    cluster_out = []
    for cluster in clusters:
        cluster_out.append({"cluster_id": cluster["cluster_id"], "period_bp": cluster["period_bp"],
                            "representative_motif": cluster["representative_motif"],
                            "array_count": cluster["array_count"], "bp": cluster["bp"],
                            "target_count": len(cluster["targets"]), "targets": ",".join(sorted(cluster["targets"])),
                            "member_count": cluster["member_count"]})
    write_tsv(out / "trf_motif_clusters.tsv", cluster_out)
    member_out = []
    for cluster in clusters:
        for member in cluster["members"]:
            member_out.append({"cluster_id": cluster["cluster_id"], "period_bp": cluster["period_bp"],
                               "representative_motif": cluster["representative_motif"], **member})
    write_tsv(out / "trf_motif_members.tsv", member_out)


def summarize_orthosite(args, out):
    rows = read_tsv(args.orthosite)
    counts = Counter((r.get("query_target", ""), r.get("evidence_class", "")) for r in rows)
    out_rows = [{"target": target, "evidence_class": evidence, "locus_count": count}
                for (target, evidence), count in sorted(counts.items())]
    write_tsv(out / "orthosite_evidence_summary.tsv", out_rows)
    family = defaultdict(lambda: Counter())
    for row in rows:
        family[(row.get("query_target", ""), row.get("pan_family", ""))][row.get("evidence_class", "")] += 1
    fam_rows = []
    for (target, pan_family), counter in sorted(family.items()):
        fam_rows.append({"target": target, "pan_family": pan_family, **dict(counter), "tested_loci": sum(counter.values())})
    write_tsv(out / "orthosite_family_evidence.tsv", fam_rows)


def parse_rmout(path: Path):
    hits = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                start, end = int(fields[5]) - 1, int(fields[6])
                if end > start:
                    hits.append((fields[4], start, end, fields[9], fields[10], int(fields[0])))
            except (ValueError, IndexError):
                continue
    return hits


def target_from_path(path: Path):
    text = str(path).lower()
    for target in TARGETS:
        if target.lower() in text:
            return target
    aliases = {"pzi_hapa1": "pziA1", "pzi_hapa2": "pziA2", "pzi_hapb1": "pziB1", "pzi_hapb2": "pziB2",
               "psp_hapa1": "pspA1", "psp_hapa2": "pspA2", "psp_hapb1": "pspB1", "psp_hapb2": "pspB2"}
    for alias, target in aliases.items():
        if alias in text:
            return target
    return ""


def summarize_nested(args, out):
    if not args.repeatmasker_root:
        return
    files = sorted(args.repeatmasker_root.rglob("*.out"))
    all_out = []
    for path in files:
        target = target_from_path(path)
        if target:
            all_out.append((target, path))
    if not all_out:
        raise FileNotFoundError(f"no target-resolvable RepeatMasker .out files under {args.repeatmasker_root}")
    rows = []
    for target, path in all_out:
        hits = parse_rmout(path)
        by_seq = defaultdict(list)
        for seqid, start, end, name, cls, score in hits:
            by_seq[seqid].append((start, end, name, cls))
        nested_bp = 0
        nested_pairs = 0
        for seqid, seq_hits in by_seq.items():
            seq_hits.sort()
            active = []
            for start, end, name, cls in seq_hits:
                active = [item for item in active if item[1] > start]
                for other_start, other_end, other_name, other_cls in active:
                    overlap = min(end, other_end) - max(start, other_start)
                    if overlap > 0 and (name != other_name or cls != other_cls):
                        nested_bp += overlap
                        nested_pairs += 1
                active.append((start, end, name, cls))
        raw_bp = sum(end - start for seqid, start, end, name, cls, score in hits)
        rows.append({"target": target, "group": GROUP[target], "source_file": path.name,
                     "raw_hit_count": len(hits), "raw_hit_bp": raw_bp,
                     "nested_overlap_bp_pair_sum": nested_bp, "nested_pair_count": nested_pairs,
                     "nested_overlap_fraction_of_raw_hits": f"{nested_bp / raw_bp:.8f}" if raw_bp else "0"})
    write_tsv(out / "te_nested_fragmentation.tsv", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drivers", type=Path, required=True)
    parser.add_argument("--divergence", type=Path, required=True)
    parser.add_argument("--intact-ltr", type=Path, required=True)
    parser.add_argument("--solo-by-target-lineage", type=Path, required=True)
    parser.add_argument("--arrays", type=Path, required=True)
    parser.add_argument("--orthosite", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeatmasker-root", type=Path)
    parser.add_argument("--species-mcmctree", type=Path)
    parser.add_argument("--subgenome-mcmctree", type=Path)
    parser.add_argument("--mu", type=float)
    parser.add_argument("--mu-time-unit-mya", type=float, default=100.0)
    parser.add_argument("--max-divergence", type=float, default=20.0)
    parser.add_argument("--top-families", type=int, default=50)
    parser.add_argument("--max-motifs", type=int, default=500)
    parser.add_argument("--motif-identity", type=float, default=0.80)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summarize_family(args, args.output_dir)
    summarize_age(args, args.output_dir)
    # Keep the validated solo-LTR source in the follow-up output directory.
    write_tsv(args.output_dir / "solo_ltr_by_target_lineage.tsv", read_tsv(args.solo_by_target_lineage))
    summarize_motifs(args, args.output_dir)
    summarize_orthosite(args, args.output_dir)
    summarize_nested(args, args.output_dir)
    manifest = {
        "mu_time_unit_mya": args.mu_time_unit_mya,
        "age_warning": "RepeatMasker divergence is a consensus-distance proxy; intact LTR divergence is preferred.",
        "orthosite_warning": "unresolved loci are not equivalent to absence.",
        "nested_warning": "nested_overlap_bp_pair_sum is a pairwise-overlap diagnostic and must not be added to genome size.",
    }
    (args.output_dir / "followup_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir), "mu": args.mu}, indent=2))


if __name__ == "__main__":
    main()
