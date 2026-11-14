#!/usr/bin/env python3
"""Summarize shared/independent POInT partition models."""

from __future__ import print_function

import argparse
import csv
from collections import Counter
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
    ap.add_argument(
        "--optimized-pillars",
        required=True,
        help="Native six-column A/B pillar matrix for exact singleton likelihoods",
    )
    ap.add_argument(
        "--process-parameters",
        type=int,
        default=2,
        help="Identifiable parameters in the fixedAB process model; excludes SwitchProb",
    )
    args = ap.parse_args()

    shared_rows = read_rows(args.shared_fit)
    candidates = [
        row
        for row in shared_rows
        if str(row.get("model", "")).startswith("reference_WGD_bias_fix")
    ]
    fixed = [row for row in candidates if "fixedAB" in str(row.get("model", ""))]
    shared_row = (fixed or candidates or [None])[0]
    if shared_row is None or not fixed or args.process_parameters != 2:
        raise SystemExit(
            "ERROR: This summary requires the fixedAB model with two process parameters; switching models are not supported"
        )
    shared_lnl = float(shared_row["log_likelihood"])

    with Path(args.optimized_pillars).open() as handle:
        native = list(csv.reader(handle, delimiter="\t"))
    if native[0] != list(SPECIES) or len(native) - 1 != args.pillars:
        raise ValueError("Pillar header/count mismatch")
    independent = {}
    independent_audit = []
    for i, species in enumerate(SPECIES):
        counts = Counter((r[i] != "NONE", r[i + 3] != "NONE") for r in native[1:])
        if counts.get((False, False), 0):
            raise ValueError("Double-absence is outside the three-state model")
        exact = sum(n * math.log(n / args.pillars) for n in counts.values())
        replicated = effective_independent_lnl(args.independent_root, species)
        if abs(exact - replicated) > 0.01:
            raise ValueError(
                "Replicated POInT singleton likelihood differs from exact multinomial"
            )
        independent[species] = exact
        independent_audit.append(
            dict(
                species=species,
                both=counts[(True, True)],
                A_only=counts[(True, False)],
                B_only=counts[(False, True)],
                exact_log_likelihood=exact,
                replicated_log_likelihood=replicated,
                identifiable_parameters=2,
            )
        )
    pair = {
        "PZI_PSP": pair_lnl(args.pair_root, "PZI_PSP"),
        "PZI_PPA": pair_lnl(args.pair_root, "PZI_PPA"),
        "PSP_PPA": pair_lnl(args.pair_root, "PSP_PPA"),
    }
    # FixedAB WGD_bias_fix has FixRate and Copy1Bias, not SwitchProb.
    # There are five branches including the root stem in a three-taxon tree,
    # three in a two-taxon tree, and one in an effective single-taxon fit.
    shared_k = args.process_parameters + 5
    pair_k = args.process_parameters + 3
    # A single terminal has three observable outcomes and exactly two free
    # probabilities. FixRate and terminal length are not separately identifiable.
    singleton_k = 2
    rows = [
        ("all_shared", "PZI+PSP+PPA shared", shared_lnl, shared_k),
        (
            "PZI_PSP_shared_PPA_independent",
            "PZI+PSP shared; PPA independent",
            pair["PZI_PSP"] + independent["PPA"],
            pair_k + singleton_k,
        ),
        (
            "PZI_PPA_shared_PSP_independent",
            "PZI+PPA shared; PSP independent",
            pair["PZI_PPA"] + independent["PSP"],
            pair_k + singleton_k,
        ),
        (
            "PSP_PPA_shared_PZI_independent",
            "PSP+PPA shared; PZI independent",
            pair["PSP_PPA"] + independent["PZI"],
            pair_k + singleton_k,
        ),
        (
            "all_independent",
            "PZI+PSP+PPA independent",
            sum(independent.values()),
            3 * singleton_k,
        ),
    ]
    n = args.pillars
    scores = []
    for code, label, lnl, k in rows:
        aic = 2.0 * k - 2.0 * lnl
        bic = k * math.log(n) - 2.0 * lnl
        scores.append(
            {
                "model": code,
                "label": label,
                "log_likelihood": lnl,
                "parameters": k,
                "pillars": n,
                "AIC": aic,
                "BIC": bic,
            }
        )
    min_aic = min(row["AIC"] for row in scores)
    min_bic = min(row["BIC"] for row in scores)
    for row in scores:
        row["delta_AIC"] = row["AIC"] - min_aic
        row["delta_BIC"] = row["BIC"] - min_bic
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with (out.parent / "singleton_identifiability_audit.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(independent_audit[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(independent_audit)
    fields = [
        "model",
        "label",
        "log_likelihood",
        "parameters",
        "pillars",
        "AIC",
        "BIC",
        "delta_AIC",
        "delta_BIC",
    ]
    with out.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in scores:
            writer.writerow(
                {
                    field: (
                        "{:.10f}".format(row[field])
                        if isinstance(row[field], float)
                        else row[field]
                    )
                    for field in fields
                }
            )
    with (out.parent / "component_likelihoods.tsv").open(
        "w", encoding="ascii", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["component", "log_likelihood"])
        writer.writerow(["shared_all", "{:.10f}".format(shared_lnl)])
        for key in ("PZI_PSP", "PZI_PPA", "PSP_PPA"):
            writer.writerow(["pair_shared_" + key, "{:.10f}".format(pair[key])])
        for species in SPECIES:
            writer.writerow(
                ["independent_" + species, "{:.10f}".format(independent[species])]
            )
    with (out.parent / "README.txt").open("w", encoding="ascii") as handle:
        handle.write(
            "Five POInT partition models fitted to the same {} optimized pillars.\n"
            "Pair models are the sum of a two-taxon shared-WGD likelihood and a one-species independent likelihood.\n"
            "Singleton likelihoods are exact multinomial log-likelihoods (three states, two identifiable parameters).\n"
            "AIC/BIC are descriptive rankings: linked pillars are not independent replicates.\n"
            "All BIC values use n={} pillars.\n".format(n, n)
        )
    print("[OK] Wrote five-model POInT comparison: {}".format(out))
    for row in scores:
        print(
            "{}\tlnL={:.6f}\tAIC={:.6f}\tBIC={:.6f}".format(
                row["model"], row["log_likelihood"], row["AIC"], row["BIC"]
            )
        )


if __name__ == "__main__":
    main()
