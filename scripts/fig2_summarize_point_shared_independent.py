#!/usr/bin/env python3
"""Compare shared-WGD and independent-polyploidy POInT fits."""

from __future__ import print_function

import argparse
import csv
import math
import re
from pathlib import Path

SPECIES = ("PZI", "PSP", "PPA")


def read_tsv(path):
    with Path(path).open("r", encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tree_ln_likelihood(path):
    text = Path(path).read_text(encoding="ascii", errors="replace")
    values = re.findall(r"Final\s+ln\s+likelihood:\s*([-+0-9.eE]+)", text, flags=re.I)
    if not values:
        raise ValueError("no Final ln likelihood in {}".format(path))
    return float(values[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared-fit", required=True, help="shared POInT point_model_fit_all.tsv")
    ap.add_argument("--independent-root", required=True)
    ap.add_argument("--pillars", required=True, type=int)
    ap.add_argument("--replicate-factor", type=int, default=3)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    shared_rows = read_tsv(args.shared_fit)
    shared = next((r for r in shared_rows if r.get("model") == "reference_WGD_bias_fix"), None)
    if shared is None:
        raise SystemExit("ERROR: shared reference_WGD_bias_fix row not found")
    shared_lnl = float(shared["log_likelihood"])
    # The r13 POInT log reports eight optimized parameters: three model
    # parameters plus five non-zero branches, including the root/stem branch
    # from the fixed Dupl state to the first species split.
    shared_k = 3 + 5

    species_rows = []
    independent_lnl = 0.0
    for species in SPECIES:
        tree = Path(args.independent_root) / species / "species_tree.nex.out"
        total = tree_ln_likelihood(tree)
        effective = total / float(args.replicate_factor)
        independent_lnl += effective
        species_rows.append((species, total, effective, str(tree)))
    # Each independent event has its own three model parameters and one
    # effective terminal branch: 3 * (3 + 1) = 12 parameters.
    independent_k = 3 * (3 + 1)
    common_n = args.pillars
    independent_n = args.pillars * len(SPECIES)

    rows = [
        {"model": "shared_WGD", "log_likelihood": shared_lnl, "parameters": shared_k,
         "sites_for_BIC": common_n, "replicate_factor": 1,
         "parameter_basis": "3 process + 5 POInT branches"},
        {"model": "independent_polyploidy", "log_likelihood": independent_lnl,
         "parameters": independent_k, "sites_for_BIC": common_n,
         "replicate_factor": args.replicate_factor,
         "parameter_basis": "3 events x (3 process + 1 terminal branch)"},
    ]
    best_aic = min(2 * r["parameters"] - 2 * r["log_likelihood"] for r in rows)
    best_bic_common = min(r["parameters"] * math.log(common_n) - 2 * r["log_likelihood"] for r in rows)
    best_bic_native = min(
        (rows[0]["parameters"] * math.log(common_n) - 2 * rows[0]["log_likelihood"]),
        (rows[1]["parameters"] * math.log(independent_n) - 2 * rows[1]["log_likelihood"]),
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as handle:
        fields = ["model", "log_likelihood", "parameters", "parameter_basis", "sites_for_BIC", "AIC", "BIC_common_n", "BIC_native_n", "delta_AIC", "delta_BIC_common_n", "delta_BIC_native_n", "replicate_factor"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in rows:
            aic = 2 * r["parameters"] - 2 * r["log_likelihood"]
            bic_common = r["parameters"] * math.log(common_n) - 2 * r["log_likelihood"]
            bic_native = r["parameters"] * math.log(common_n if r["model"] == "shared_WGD" else independent_n) - 2 * r["log_likelihood"]
            writer.writerow({
                "model": r["model"], "log_likelihood": "{:.10f}".format(r["log_likelihood"]),
                "parameters": r["parameters"], "parameter_basis": r["parameter_basis"], "sites_for_BIC": r["sites_for_BIC"],
                "AIC": "{:.6f}".format(aic), "BIC_common_n": "{:.6f}".format(bic_common),
                "BIC_native_n": "{:.6f}".format(bic_native),
                "delta_AIC": "{:.6f}".format(aic - best_aic),
                "delta_BIC_common_n": "{:.6f}".format(bic_common - best_bic_common),
                "delta_BIC_native_n": "{:.6f}".format(bic_native - best_bic_native),
                "replicate_factor": r["replicate_factor"],
            })
    species_out = out.with_name("independent_species_likelihoods.tsv")
    with species_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "replicated_total_lnL", "effective_lnL", "tree_output"])
        writer.writerows(species_rows)
    print("[OK] Shared/independent POInT comparison written: {}".format(out))
    print("[INFO] shared lnL={:.6f}; independent summed effective lnL={:.6f}".format(shared_lnl, independent_lnl))


if __name__ == "__main__":
    main()
