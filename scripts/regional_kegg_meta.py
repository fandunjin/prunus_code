#!/usr/bin/env python3
"""Cross-species meta-analysis and pathway-class enrichment for non-overlapping A/B regions.

Mirrors the GO workflow (regional_go_meta.py): combines per-species
hypergeometric p-values by a Bonferroni partial-conjunction test (at least two of three species) for the same KEGG pathway, and tests
KEGG BRITE class B (pathway category) enrichment per species, then across
species, so shared enrichment across PZI, PSP and PPA can be assessed directly.
"""

from __future__ import annotations

from enrichment_utils import bh_group, hypergeom_sf, partial_conjunction

from collections import defaultdict
from pathlib import Path

import pandas as pd


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
ROOT = args.root
OUT = ROOT / "GO"
ANN = ROOT / "annotation" / "functional"
SPECIES = ["PZI", "PSP", "PPA"]
SHORT = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}


def read_pathways(species: str) -> pd.DataFrame:
    return pd.read_csv(
        ANN / species / "swissprot_kegg_pathway.tsv",
        sep="\t",
        dtype=str,
        usecols=[
            "protein_id",
            "pathway_id",
            "pathway_name",
            "pathway_class_a",
            "pathway_class_b",
        ],
    ).dropna(subset=["protein_id", "pathway_id"])


def main() -> None:
    members = pd.read_csv(OUT / "KEGG_nonoverlap_gene_members.tsv", sep="\t")
    ora = pd.read_csv(OUT / "KEGG_nonoverlap_ORA.tsv", sep="\t")
    pathways_by_species = {sp: read_pathways(SHORT[sp]) for sp in SPECIES}

    meta_rows = []
    theme_rows = []
    for direction in ("A_favoured", "B_favoured"):
        sub = ora[
            ora["scope"].eq("shared_only")
            & ora["background"].eq("genome")
            & ora["direction"].eq(direction)
            & ora["term_background"].ge(10)
        ]
        for (pathway_id, name, class_a, class_b), grp in sub.groupby(
            ["pathway_id", "pathway_name", "pathway_class_a", "pathway_class_b"]
        ):
            if grp["species"].nunique() < 2:
                continue
            meta_rows.append(
                {
                    "direction": direction,
                    "pathway_id": pathway_id,
                    "pathway_name": name,
                    "pathway_class_a": class_a,
                    "pathway_class_b": class_b,
                    "n_species_tested": grp["species"].nunique(),
                    "n_species_fdr05": int(grp["fdr"].lt(0.05).sum()),
                    "species": ",".join(sorted(grp["species"])),
                    "mean_fold": float(grp["fold_enrichment"].mean()),
                    "min_fdr": float(grp["fdr"].min()),
                    "combined_p": partial_conjunction(grp["p_value"].tolist()),
                }
            )

        fg_genes = {
            species: set(
                members[
                    members["species"].eq(species) & members["direction"].eq(direction)
                ]["gene_id"]
                .str.split("__", n=1)
                .str[-1]
            )
            for species in SPECIES
        }
        for species in SPECIES:
            table = pathways_by_species[species]
            all_genes = set(table["protein_id"])
            fg = fg_genes[species] & all_genes
            if not fg:
                continue
            class_genes: dict[str, set] = defaultdict(set)
            for protein_id, class_b in table.dropna(subset=["pathway_class_b"])[
                ["protein_id", "pathway_class_b"]
            ].itertuples(index=False):
                if class_b:
                    class_genes[class_b].add(protein_id)
            n, N = len(fg), len(all_genes)
            for theme, theme_genes in class_genes.items():
                k = len(fg & theme_genes)
                K = len(theme_genes)
                theme_rows.append(
                    {
                        "direction": direction,
                        "species": species,
                        "theme": theme,
                        "foreground_genes": n,
                        "background_genes": N,
                        "term_foreground": k,
                        "term_background": K,
                        "fold_enrichment": (k / n) / (K / N),
                        "p_value": hypergeom_sf(k, N, K, n),
                    }
                )

    meta = pd.DataFrame(meta_rows)
    if not meta.empty:
        meta["fdr"] = meta.groupby(["direction"])["combined_p"].transform(bh_group)
        meta.sort_values(["direction", "combined_p"]).to_csv(
            OUT / "KEGG_nonoverlap_meta.tsv", sep="\t", index=False
        )

    themes = pd.DataFrame(theme_rows)
    if not themes.empty:
        themes["fdr"] = themes.groupby(["direction", "theme"])["p_value"].transform(
            bh_group
        )
        meta_theme = (
            themes.groupby(["direction", "theme"])
            .agg(
                n_species_tested=("species", "nunique"),
                n_species_fdr05=("fdr", lambda v: int((v < 0.05).sum())),
                mean_fold=("fold_enrichment", "mean"),
                min_fdr=("fdr", "min"),
                combined_p=("p_value", lambda v: partial_conjunction(list(v))),
            )
            .reset_index()
        )
        meta_theme["fdr_combined"] = meta_theme.groupby(["direction"])[
            "combined_p"
        ].transform(bh_group)
        themes.to_csv(OUT / "KEGG_nonoverlap_theme.tsv", sep="\t", index=False)
        meta_theme.to_csv(OUT / "KEGG_nonoverlap_theme_meta.tsv", sep="\t", index=False)

    if not meta.empty:
        print("== pathway-level meta ==")
        print(
            meta.sort_values("combined_p")
            .head(15)[
                [
                    "direction",
                    "pathway_name",
                    "pathway_class_b",
                    "n_species_tested",
                    "n_species_fdr05",
                    "mean_fold",
                    "combined_p",
                    "fdr",
                ]
            ]
            .to_string(index=False)
        )
    if "meta_theme" in dir():
        print()
        print("== pathway-class (theme) level ==")
        print(meta_theme.sort_values("combined_p").to_string(index=False))


if __name__ == "__main__":
    main()
