#!/usr/bin/env python3
"""Reconstruct candidate anchor attrition with Dollo parsimony for Fig. 2d.

The input is the conservative GTR-centered JCVI state table.  Loci are retained
when PZI, PSP and PPA each have at least one unambiguous A/B anchor.  Starting with both copies after the
registered duplication, losses are placed on the smallest compatible branch.
The output is anchor attrition, not POInT likelihood-based gene loss.
"""

from __future__ import print_function

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


SPECIES = ("PZI", "PSP", "PPA")
SUBGENOMES = ("A", "B")
BRANCHES = ("post_duplication", "prunus_stem", "psp_pzi_stem", "pzi_tip", "psp_tip", "ppa_tip")


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def state_to_presence(state):
    if state == "both_anchored":
        return {"A": 1, "B": 1}
    if state == "A_only_anchored":
        return {"A": 1, "B": 0}
    if state == "B_only_anchored":
        return {"A": 0, "B": 1}
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    by_gene = defaultdict(dict)
    with Path(args.states).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("dataset") != "main" or row.get("species") not in SPECIES:
                continue
            gene = row.get("reference_gene", "")
            if gene:
                by_gene[gene][row["species"]] = row
    eligible = []
    excluded = []
    for gene, species_rows in sorted(by_gene.items()):
        if set(species_rows) != set(SPECIES):
            excluded.append({"reference_gene": gene, "reason": "missing_species_state"})
            continue
        presence = {species: state_to_presence(species_rows[species]["anchor_state"]) for species in SPECIES}
        if any(value is None for value in presence.values()):
            excluded.append({"reference_gene": gene, "reason": "ambiguous_or_neither_anchored"})
            continue
        eligible.append((gene, species_rows, presence))
    if not eligible:
        raise SystemExit("ERROR: no eligible conservative pillars for JCVI-Dollo reconstruction")

    branch_states = defaultdict(dict)
    branch_losses = Counter()
    locus_rows = []
    for gene, species_rows, presence in eligible:
        for subgenome in SUBGENOMES:
            tips = {species: presence[species][subgenome] for species in SPECIES}
            post = 1
            prunus = 1 if any(tips.values()) else 0
            psp_pzi = 1 if tips["PZI"] or tips["PSP"] else 0
            states = {
                "post_duplication": post,
                "prunus_stem": prunus,
                "psp_pzi_stem": psp_pzi,
                "pzi_tip": tips["PZI"],
                "psp_tip": tips["PSP"],
                "ppa_tip": tips["PPA"],
            }
            losses = {
                "post_duplication": 0,
                "prunus_stem": post - prunus,
                "psp_pzi_stem": max(0, prunus - psp_pzi),
                "pzi_tip": max(0, psp_pzi - tips["PZI"]),
                "psp_tip": max(0, psp_pzi - tips["PSP"]),
                "ppa_tip": max(0, prunus - tips["PPA"]),
            }
            for branch in BRANCHES:
                branch_states[(branch, subgenome)][gene] = states[branch]
                branch_losses[(branch, subgenome)] += losses[branch]
            locus_rows.append({
                "reference_gene": gene,
                "reference_seqid": species_rows["PZI"].get("reference_seqid", ""),
                "subgenome": subgenome,
                "pzi_state": tips["PZI"],
                "psp_state": tips["PSP"],
                "ppa_state": tips["PPA"],
                "prunus_inferred_state": prunus,
                "psp_pzi_inferred_state": psp_pzi,
            })

    parent = {
        "post_duplication": None,
        "prunus_stem": "post_duplication",
        "psp_pzi_stem": "prunus_stem",
        "pzi_tip": "psp_pzi_stem",
        "psp_tip": "psp_pzi_stem",
        "ppa_tip": "prunus_stem",
    }
    summary_rows = []
    for branch in BRANCHES:
        for subgenome in SUBGENOMES:
            surviving = sum(branch_states[(branch, subgenome)].values())
            parent_branch = parent[branch]
            parent_surviving = len(eligible) if parent_branch is None else sum(branch_states[(parent_branch, subgenome)].values())
            ratio = float(surviving) / parent_surviving if parent_surviving else 0.0
            summary_rows.append({
                "branch_id": branch,
                "subgenome": subgenome,
                "surviving_anchors": surviving,
                "anchor_attrition": branch_losses[(branch, subgenome)],
                "retention_ratio": "{:.8f}".format(ratio),
                "n_eligible_pillars": len(eligible),
                "analysis_status": "complete",
                "method": "JCVI_anchored_Dollo_parsimony",
            })
    write_tsv(outdir / "jcvi_dollo_branch_summary.tsv", list(summary_rows[0].keys()), summary_rows)
    write_tsv(outdir / "jcvi_dollo_locus_states.tsv", list(locus_rows[0].keys()), locus_rows)
    write_tsv(outdir / "excluded_loci.tsv", ["reference_gene", "reason"], excluded)
    (outdir / "method_status.txt").write_text(
        "Dollo parsimony was applied to conservative JCVI anchor presence on (GTR,(PPA,(PSP,PZI))).\n"
        "Numbers are surviving syntenic anchors and candidate anchor attrition, not POInT gene-loss likelihood estimates.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("JCVI-anchored Dollo reconstruction completed.\n", encoding="ascii")
    print("[OK] JCVI-Dollo eligible loci={}; excluded loci={}".format(len(eligible), len(excluded)))


if __name__ == "__main__":
    main()
