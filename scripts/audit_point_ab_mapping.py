#!/usr/bin/env python3
"""Audit physical A/B states against POInT latent tracking states.

POInT's Copy1 and Copy2 labels are latent states.  The posterior tracking
pattern (for example PZI01 or PZI10) records how the physical A/B inputs map
to those states at each pillar.  This script keeps those quantities separate
and writes a compact, source-data-traceable comparison.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SPECIES = ("PZI", "PSP", "PPA")
PATTERN_RE = re.compile(r"^(PZI|PSP|PPA)([01]{2})$")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def present(value):
    return bool((value or "").strip()) and (value or "").strip().upper() != "NONE"


def pct(value, total):
    return 100.0 * value / total if total else 0.0


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pillars", required=True, help="GTR-centered pillar TSV")
    parser.add_argument("--posterior", required=True, help="normalized POInT posterior summary TSV")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    pillars = read_tsv(args.pillars)
    posterior = read_tsv(args.posterior)
    if len(pillars) != len(posterior):
        raise SystemExit(
            "pillar/posterior count mismatch: {} versus {}".format(len(pillars), len(posterior))
        )
    if not pillars:
        raise SystemExit("empty pillar table")

    # POInT writes posterior genes in its current tracking-slot order, which
    # can swap relative to physical A/B.  Verify row alignment by comparing
    # the unordered non-missing gene IDs rather than slot positions.
    posterior_fields = [f"{sp}_{copy}_gene" for sp in SPECIES for copy in ("A", "B")]
    for index, (pillar, post) in enumerate(zip(pillars, posterior)):
        expected_genes = {
            value
            for field in posterior_fields
            for value in [pillar.get(field, "") or ""]
            if present(value)
        }
        observed_genes = {
            value
            for field in posterior_fields
            for value in [post.get(field, "") or ""]
            if present(value)
        }
        if expected_genes != observed_genes:
            raise SystemExit(
                "gene-set mismatch at pillar {}: {} != {}".format(
                    index, sorted(expected_genes), sorted(observed_genes)
                )
            )

    rows = []
    mapping_rows = []
    for species in SPECIES:
        direct_a = sum(present(row.get(f"{species}_A_gene")) for row in pillars)
        direct_b = sum(present(row.get(f"{species}_B_gene")) for row in pillars)
        both = sum(
            present(row.get(f"{species}_A_gene")) and present(row.get(f"{species}_B_gene"))
            for row in pillars
        )
        a_only = sum(
            present(row.get(f"{species}_A_gene")) and not present(row.get(f"{species}_B_gene"))
            for row in pillars
        )
        b_only = sum(
            present(row.get(f"{species}_B_gene")) and not present(row.get(f"{species}_A_gene"))
            for row in pillars
        )
        pattern_counts = {"01": 0, "10": 0}
        latent_c1 = 0
        latent_c2 = 0
        for pillar, post in zip(pillars, posterior):
            token = post.get("best_tracking", "").split("_")
            token = next((part for part in token if part.startswith(species)), "")
            match = PATTERN_RE.fullmatch(token)
            if not match:
                raise SystemExit("invalid tracking state at pillar {}: {}".format(len(mapping_rows), token))
            pattern = match.group(2)
            pattern_counts[pattern] += 1
            physical_a = present(pillar.get(f"{species}_A_gene"))
            physical_b = present(pillar.get(f"{species}_B_gene"))
            # 01 means physical A -> latent Copy1 and physical B -> Copy2;
            # 10 swaps those assignments for this pillar.
            if pattern == "01":
                latent_c1 += int(physical_a)
                latent_c2 += int(physical_b)
            else:
                latent_c1 += int(physical_b)
                latent_c2 += int(physical_a)
        rows.extend(
            [
                {
                    "species": species,
                    "quantity": "physical_A_retained",
                    "count": direct_a,
                    "percent_of_13204": "{:.4f}".format(pct(direct_a, len(pillars))),
                    "definition": "observed physical A input presence",
                },
                {
                    "species": species,
                    "quantity": "physical_B_retained",
                    "count": direct_b,
                    "percent_of_13204": "{:.4f}".format(pct(direct_b, len(pillars))),
                    "definition": "observed physical B input presence",
                },
                {
                    "species": species,
                    "quantity": "both_physical_copies",
                    "count": both,
                    "percent_of_13204": "{:.4f}".format(pct(both, len(pillars))),
                    "definition": "both physical A and B present",
                },
                {
                    "species": species,
                    "quantity": "A_only_physical",
                    "count": a_only,
                    "percent_of_13204": "{:.4f}".format(pct(a_only, len(pillars))),
                    "definition": "physical A present and physical B absent",
                },
                {
                    "species": species,
                    "quantity": "B_only_physical",
                    "count": b_only,
                    "percent_of_13204": "{:.4f}".format(pct(b_only, len(pillars))),
                    "definition": "physical B present and physical A absent",
                },
                {
                    "species": species,
                    "quantity": "latent_Copy1_best_posterior_presence",
                    "count": latent_c1,
                    "percent_of_13204": "{:.4f}".format(pct(latent_c1, len(pillars))),
                    "definition": "presence after applying each pillar's best POInT tracking permutation",
                },
                {
                    "species": species,
                    "quantity": "latent_Copy2_best_posterior_presence",
                    "count": latent_c2,
                    "percent_of_13204": "{:.4f}".format(pct(latent_c2, len(pillars))),
                    "definition": "presence after applying each pillar's best POInT tracking permutation",
                },
                {
                    "species": species,
                    "quantity": "posterior_pattern_01",
                    "count": pattern_counts["01"],
                    "percent_of_13204": "{:.4f}".format(pct(pattern_counts["01"], len(pillars))),
                    "definition": "physical A maps to Copy1 and physical B maps to Copy2",
                },
                {
                    "species": species,
                    "quantity": "posterior_pattern_10",
                    "count": pattern_counts["10"],
                    "percent_of_13204": "{:.4f}".format(pct(pattern_counts["10"], len(pillars))),
                    "definition": "physical B maps to Copy1 and physical A maps to Copy2",
                },
            ]
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fields = ["species", "quantity", "count", "percent_of_13204", "definition"]
    write_tsv(outdir / "point_ab_mapping_audit.tsv", rows, fields)
    (outdir / "README.txt").write_text(
        "The 13,204 pillars were checked for identical physical A/B genes before comparison.\n"
        "Physical A/B retention is an observed anchor-presence summary. Copy1/Copy2 are POInT latent\n"
        "states; their correspondence to physical A/B is pillar-specific and comes from best_tracking.\n"
        "Best-posterior latent presence is descriptive and is not a substitute for POInT model-expected\n"
        "retention or an independently validated gene-loss count.\n",
        encoding="ascii",
    )
    print("[OK] Wrote POInT physical/latent copy audit to {}".format(outdir))


if __name__ == "__main__":
    main()
