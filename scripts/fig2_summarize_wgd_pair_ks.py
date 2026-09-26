#!/usr/bin/env python3
"""Audit and filter pairwise WGD v2 CodeML Ks output for Figure 2a."""

from __future__ import print_function

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ks", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--max-ks", type=float, default=5.0)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    families = pd.read_csv(args.families, sep="\t", dtype=str)
    if "family_id" not in families.columns:
        raise SystemExit("[ERROR] Missing family_id in {}".format(args.families))
    gene_fields = [column for column in families.columns if column != "family_id"]
    if len(gene_fields) != 2:
        raise SystemExit("[ERROR] Expected two gene columns; found {}".format(gene_fields))
    if families["family_id"].duplicated().any():
        raise SystemExit("[ERROR] Duplicate family_id values in {}".format(args.families))

    ks = pd.read_csv(args.ks, sep="\t", dtype={"family": str})
    required = {"family", "gene1", "gene2", "dS"}
    missing = required - set(ks.columns)
    if missing:
        raise SystemExit("[ERROR] Missing WGD columns: {}".format(",".join(sorted(missing))))
    if ks["family"].duplicated().any():
        duplicates = ks.loc[ks["family"].duplicated(), "family"].head().tolist()
        raise SystemExit("[ERROR] More than one Ks row per pair family: {}".format(duplicates))
    unexpected = sorted(set(ks["family"].dropna()) - set(families["family_id"]))
    if unexpected:
        raise SystemExit("[ERROR] Unexpected WGD families: {}".format(unexpected[:5]))

    merged = families.merge(ks, how="left", left_on="family_id", right_on="family", validate="one_to_one")
    present = merged["family"].notna()
    gene_match = []
    for _, row in merged.iterrows():
        if pd.isna(row.get("family")):
            gene_match.append(False)
            continue
        expected = {str(row[gene_fields[0]]), str(row[gene_fields[1]])}
        observed = {str(row["gene1"]), str(row["gene2"])}
        gene_match.append(expected == observed)
    merged["gene_ids_match_family"] = gene_match
    mismatched = present & ~merged["gene_ids_match_family"]
    if mismatched.any():
        bad = merged.loc[mismatched, "family_id"].head().tolist()
        raise SystemExit("[ERROR] WGD gene IDs do not match input pair families: {}".format(bad))

    merged["dS"] = pd.to_numeric(merged["dS"], errors="coerce")
    finite = np.isfinite(merged["dS"])
    merged["ks_status"] = "retained_0_lt_Ks_lt_{}".format(args.max_ks)
    merged.loc[~present, "ks_status"] = "missing_wgd_result"
    merged.loc[present & ~finite, "ks_status"] = "nonfinite_Ks"
    merged.loc[present & finite & (merged["dS"] <= 0), "ks_status"] = "nonpositive_Ks"
    merged.loc[present & finite & (merged["dS"] >= args.max_ks), "ks_status"] = "saturated_Ks_ge_{}".format(args.max_ks)
    retained_mask = present & finite & (merged["dS"] > 0) & (merged["dS"] < args.max_ks)

    merged.insert(0, "dataset", args.dataset)
    retained = merged.loc[retained_mask].copy()
    excluded = merged.loc[~retained_mask].copy()
    merged.to_csv(args.outdir / (args.dataset + ".all_pair_ks_audit.tsv"), sep="\t", index=False, na_rep="NA")
    retained.to_csv(args.outdir / (args.dataset + ".ks_0_5.tsv"), sep="\t", index=False, na_rep="NA")
    excluded.to_csv(args.outdir / (args.dataset + ".ks_excluded.tsv"), sep="\t", index=False, na_rep="NA")

    status_counts = merged["ks_status"].value_counts().to_dict()
    summary = pd.DataFrame([{
        "dataset": args.dataset,
        "families_input": len(families),
        "wgd_rows": int(present.sum()),
        "retained_0_lt_Ks_lt_5": int(retained_mask.sum()),
        "missing_wgd_result": int((merged["ks_status"] == "missing_wgd_result").sum()),
        "nonfinite_Ks": int((merged["ks_status"] == "nonfinite_Ks").sum()),
        "nonpositive_Ks": int((merged["ks_status"] == "nonpositive_Ks").sum()),
        "saturated_Ks_ge_5": int((merged["ks_status"] == "saturated_Ks_ge_5.0").sum()),
        "retention_percent": 100.0 * retained_mask.sum() / len(families),
        "Ks_min": retained["dS"].min(),
        "Ks_median": retained["dS"].median(),
        "Ks_mean": retained["dS"].mean(),
        "Ks_max": retained["dS"].max(),
    }])
    summary.to_csv(args.outdir / (args.dataset + ".ks_summary.tsv"), sep="\t", index=False, na_rep="NA")
    (args.outdir / "COMPLETE").write_text("WGD pairwise Ks audit completed.\n", encoding="ascii")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
