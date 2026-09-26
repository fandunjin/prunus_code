#!/usr/bin/env python3
"""Integrate neutral internal synteny evidence into PZI chromosome models.

This is an interpretation guardrail, not a new likelihood calculation.  It
prevents GTR multi-block anchors from being counted as fusion/fission evidence
when the same copy1 chromosome has an independently confirmed primary copy2
homolog over most of its sequence.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model-comparison", required=True)
    p.add_argument("--k97-partitions", required=True)
    p.add_argument("--resolution", required=True)
    p.add_argument("--hap-pairs", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read(path):
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, fields, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def main():
    a = args(); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    models = read(a.model_comparison)
    resolution = read(a.resolution)
    hap_pairs = read(a.hap_pairs)
    top_k97 = read(a.k97_partitions)
    if not resolution:
        raise SystemExit("No copy1/copy2 resolution rows")
    n = len(resolution)
    primary_top = sum(int(r["expected_target_is_top"]) for r in resolution if r["expected_target_is_top"] != "NA")
    expected_fractions = [float(r["expected_target_anchor_fraction"]) for r in resolution if r["expected_target_anchor_fraction"] != "NA"]
    genome_fractions = [float(r["same_chr_genome_alignment_fraction"]) for r in resolution]
    gtr_multi = [r for r in resolution if "," in r["gtr_supported_blocks"]]
    multi_primary = sum(r["interpretation"] == "GTR_multiblock_stable_primary_internal_homolog" for r in gtr_multi)
    candidates = [r for r in resolution if r["interpretation"] == "candidate_copy_specific_breakpoint"]
    hap_ranked = sorted(hap_pairs, key=lambda r: (float(r["query_top_specificity"]) + float(r["reference_top_specificity"])) / 2, reverse=True)
    top_hap = hap_ranked[:2]
    top_text = "; ".join(
        "{} (mean chromosome specificity {:.3f}; 8x8)".format(
            r["comparison"], (float(r["query_top_specificity"]) + float(r["reference_top_specificity"])) / 2)
        for r in top_hap)
    partition_nonunique = "NA"
    if len(top_k97) >= 2:
        partition_nonunique = int(top_k97[0]["score"] == top_k97[1]["score"])

    evidence = [{
        "metric": "copy1_chromosomes_tested", "value": n,
        "interpretation": "All normalized copy1 chr01-chr16 were tested against their normalized copy2 counterpart."
    }, {
        "metric": "expected_copy2_target_is_top_anchor_target", "value": f"{primary_top}/{n}",
        "interpretation": "Gene-level primary targets agree with the normalized copy1/copy2 homolog map."
    }, {
        "metric": "expected_copy2_anchor_fraction_mean", "value": f"{statistics.mean(expected_fractions):.6f}",
        "interpretation": "Remaining anchors are expected homeologous/duplicated matches and cannot alone indicate fission."
    }, {
        "metric": "same_chromosome_genome_alignment_fraction_median", "value": f"{statistics.median(genome_fractions):.6f}",
        "interpretation": "Genome alignment independently validates the primary same-chromosome relation."
    }, {
        "metric": "GTR_multiblock_copy1_chromosomes", "value": len(gtr_multi),
        "interpretation": "These chromosomes carry anchors to more than one GTR block."
    }, {
        "metric": "GTR_multiblock_with_primary_internal_homolog", "value": f"{multi_primary}/{len(gtr_multi)}",
        "interpretation": "GTR multi-block anchors do not independently establish copy-specific fusions/fissions."
    }, {
        "metric": "copy_specific_breakpoint_candidates", "value": len(candidates),
        "interpretation": "Candidates require discordant primary copy2 target or weak same-chromosome genome alignment."
    }, {
        "metric": "two_highest_independent_haplotype_pairs", "value": top_text,
        "interpretation": "Labels are neutral; these are the two data-driven strongest 8x8 pairings among four haplotypes."
    }, {
        "metric": "K9plus7_top_partition_tie", "value": partition_nonunique,
        "interpretation": "A tie means the 9/7 partition is not uniquely identified by the reference-anchor score."
    }]
    write(out / "internal_evidence_summary.tsv", ["metric", "value", "interpretation"], evidence)

    by_model = {r["model"]: r for r in models}
    rows = [
        {"model": "K=8", "reference_score": by_model.get("K=8", {}).get("score", "NA"),
         "neutral_internal_evidence": "Two data-driven strongest haplotype pairings are both 8x8; 16/16 copy1 chromosomes recover their primary copy2 homolog.",
         "model_status": "consistent but not alone decisive",
         "claim_boundary": "Do not treat the original GTR-profile K=8 pair list as final homolog pairs without direct internal block validation."},
        {"model": "K=9-block", "reference_score": by_model.get("K=9-block", {}).get("score", "NA"),
         "neutral_internal_evidence": f"{multi_primary}/{len(gtr_multi)} GTR-multiblock chromosomes retain a primary same-chromosome copy2 homolog; no copy-specific breakpoint candidate was recovered.",
         "model_status": "GTR block structure remains possible; fusion/fission count is not independently supported",
         "claim_boundary": "Do not call the reference multi-block events true PZI fusions/fissions without continuous breakpoint and orthogonal assembly evidence."},
        {"model": "K=9+7", "reference_score": by_model.get("K=9+7", {}).get("score", "NA"),
         "neutral_internal_evidence": "The leading 9/7 reference-anchor partitions are tied/non-unique; neutral independent haplotype comparisons identify two strongest 8x8 pairings.",
         "model_status": "not independently supported",
         "claim_boundary": "A 9+7 chromosome-number assertion requires an independently resolved 9/7 partition and breakpoint evidence."},
    ]
    write(out / "model_assessment_after_internal_synteny.tsv",
          ["model", "reference_score", "neutral_internal_evidence", "model_status", "claim_boundary"], rows)
    (out / "README.txt").write_text(
        "Internal-evidence integration of neutral chromosome models. This does not replace\n"
        "model likelihoods. It restricts interpretation of GTR multi-block anchors after\n"
        "copy1-copy2 genome and gene-level homolog validation.\n", encoding="ascii")
    print(f"[OK] Internal evidence integrated: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
