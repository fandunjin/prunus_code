#!/usr/bin/env python3
"""Summarize shared/independent POInT partition models."""

from __future__ import print_function

import argparse
import csv
import math
import re
from pathlib import Path


MODELS = (
    ("all_shared", "PZI+PSP+PPA shared"),
    ("PZI_PSP_shared_PPA_independent", "PZI+PSP shared; PPA independent"),
    ("PZI_PPA_shared_PSP_independent", "PZI+PPA shared; PSP independent"),
    ("PSP_PPA_shared_PZI_independent", "PSP+PPA shared; PZI independent"),
    ("all_independent", "PZI+PSP+PPA independent"),
)
SPECIES = ("PZI", "PSP", "PPA")


def final_lnl(path):
    text = Path(path).read_text(encoding="ascii", errors="replace")
    values = re.findall(r"Final\s+ln\s+likelihood:\s*([-+0-9.eE]+)", text, flags=re.I)
    if not values:
        raise ValueError("Final ln likelihood not found in {}".format(path))
    return float(values[-1])


def effective_independent_lnl(independent_root, species):
    return final_lnl(Path(independent_root) / species / "species_tree.nex.out") / 3.0


def read_rows(path):
    with Path(path).open("r", encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pair_lnl(pair_root, pair):
    return final_lnl(Path(pair_root) / pair / "species_tree.nex.out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared-fit", required=True)
    ap.add_argument("--pair-root", required=True)
    ap.add_argument("--independent-root", required=True)
    ap.add_argument("--pillars", required=True, type=int)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    shared_rows = read_rows(args.shared_fit)
    candidates = [row for row in shared_rows
                  if str(row.get("model", "")).startswith("reference_WGD_bias_fix")]
    fixed = [row for row in candidates if "fixedAB" in str(row.get("model", ""))]
    shared_row = (fixed or candidates or [None])[0]
    if shared_row is None:
        raise SystemExit("ERROR: reference_WGD_bias_fix shared fit not found")
    shared_lnl = float(shared_row["log_likelihood"])

    independent = {species: effective_independent_lnl(args.independent_root, species) for species in SPECIES}
    pair = {
        "PZI_PSP": pair_lnl(args.pair_root, "PZI_PSP"),
        "PZI_PPA": pair_lnl(args.pair_root, "PZI_PPA"),
        "PSP_PPA": pair_lnl(args.pair_root, "PSP_PPA"),
    }
    # POInT parameters: all-shared = 3 process + 5 rooted three-taxon branches;
    # a two-taxon shared fit = 3 process + 3 branches; each independent event
    # contributes 3 process parameters + 1 effective terminal branch.
    rows = [
        ("all_shared", "PZI+PSP+PPA shared", shared_lnl, 8),
        ("PZI_PSP_shared_PPA_independent", "PZI+PSP shared; PPA independent", pair["PZI_PSP"] + independent["PPA"], 10),
        ("PZI_PPA_shared_PSP_independent", "PZI+PPA shared; PSP independent", pair["PZI_PPA"] + independent["PSP"], 10),
        ("PSP_PPA_shared_PZI_independent", "PSP+PPA shared; PZI independent", pair["PSP_PPA"] + independent["PZI"], 10),
        ("all_independent", "PZI+PSP+PPA independent", sum(independent.values()), 12),
    ]
    n = args.pillars
    scores = []
    for code, label, lnl, k in rows:
        aic = 2.0 * k - 2.0 * lnl
        bic = k * math.log(n) - 2.0 * lnl
        scores.append({"model": code, "label": label, "log_likelihood": lnl,
                       "parameters": k, "pillars": n, "AIC": aic, "BIC": bic})
    min_aic = min(row["AIC"] for row in scores)
    min_bic = min(row["BIC"] for row in scores)
    for row in scores:
        row["delta_AIC"] = row["AIC"] - min_aic
        row["delta_BIC"] = row["BIC"] - min_bic
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["model", "label", "log_likelihood", "parameters", "pillars", "AIC", "BIC", "delta_AIC", "delta_BIC"]
    with out.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in scores:
            writer.writerow({field: ("{:.10f}".format(row[field]) if isinstance(row[field], float) else row[field]) for field in fields})
    with (out.parent / "component_likelihoods.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["component", "log_likelihood"])
        writer.writerow(["shared_all", "{:.10f}".format(shared_lnl)])
        for key in ("PZI_PSP", "PZI_PPA", "PSP_PPA"):
            writer.writerow(["pair_shared_" + key, "{:.10f}".format(pair[key])])
        for species in SPECIES:
            writer.writerow(["independent_" + species, "{:.10f}".format(independent[species])])
    with (out.parent / "README.txt").open("w", encoding="ascii") as handle:
        handle.write(
            "Five POInT partition models fitted to the same {} optimized pillars.\n"
            "Pair models are the sum of a two-taxon shared-WGD likelihood and a one-species independent likelihood.\n"
            "All BIC values use n={} pillars.\n".format(n, n)
        )
    print("[OK] Wrote five-model POInT comparison: {}".format(out))
    for row in scores:
        print("{}\tlnL={:.6f}\tAIC={:.6f}\tBIC={:.6f}".format(row["model"], row["log_likelihood"], row["AIC"], row["BIC"]))


if __name__ == "__main__":
    main()
