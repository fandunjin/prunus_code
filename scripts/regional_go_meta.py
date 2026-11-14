#!/usr/bin/env python3
"""Cross-species meta-analysis and theme-level enrichment for non-overlapping A/B regions.

Combines per-species hypergeometric p-values by a Bonferroni partial-conjunction test (at least two of three species) (same GO term or
same functional theme) and reports species support, so shared enrichment across
PZI, PSP and PPA can be assessed directly.
"""

from __future__ import annotations

from enrichment_utils import bh_group, hypergeom_sf, partial_conjunction

from pathlib import Path

import pandas as pd


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
ROOT = args.root
OUT = ROOT / "GO"
GO_DIR = ROOT / "annotation" / "functional"
SPECIES = ["PZI", "PSP", "PPA"]
SHORT = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
EVIDENCE_MODES = ("all", "non_IEA")

THEMES = {
    "Transcriptional regulation": (
        "transcription",
        "rna polymerase ii",
        "cis-regulatory",
        "nucleus",
        "dna-binding transcription",
    ),
    "Cell wall & pectin": ("cell wall", "polygalacturonase", "pectin", "cuticle"),
    "Senescence": ("senescence",),
    "Development & pollen": ("pollen", "embryonic", "pattern formation"),
    "Transport": ("transporter", "antiporter", "efflux", "transmembrane transport"),
    "Hormone & ripening": ("gibberellin", "ripening", "auxin", "ethylene", "abscisic"),
    "Endomembrane & vacuole": ("endosome", "vacuole", "vacuolar"),
}


def read_go(path: Path, mode: str) -> pd.DataFrame:
    usecols = ["protein_id", "go_id", "go_aspect", "go_description", "go_evidence"]
    go = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str)
    go = go.dropna(subset=["protein_id", "go_id", "go_aspect"])
    if mode == "non_IEA":
        go = go[go["go_evidence"].fillna("").str.upper().ne("IEA")]
    return go.drop_duplicates(["protein_id", "go_id"])


def main() -> None:
    members = pd.read_csv(OUT / "GO_nonoverlap_gene_members.tsv", sep="\t")

    meta_rows = []
    theme_rows = []
    for mode in EVIDENCE_MODES:
        go_by_species = {}
        for species in SPECIES:
            go = read_go(GO_DIR / SHORT[species] / "swissprot_go.tsv", mode)
            go_by_species[species] = go

        for direction in ("A_favoured", "B_favoured"):
            # term-level meta on the per-species ORA results
            ora = pd.read_csv(OUT / "GO_nonoverlap_ORA.tsv", sep="\t")
            sub = ora[
                ora["scope"].eq("shared_only")
                & ora["background"].eq("genome")
                & ora["evidence_mode"].eq(mode)
                & ora["direction"].eq(direction)
                & ora["term_background"].ge(10)
            ]
            for (go_id, aspect, desc), grp in sub.groupby(
                ["go_id", "go_aspect", "go_description"]
            ):
                if grp["species"].nunique() < 2:
                    continue
                meta_rows.append(
                    {
                        "evidence_mode": mode,
                        "direction": direction,
                        "go_id": go_id,
                        "go_aspect": aspect,
                        "go_description": desc,
                        "n_species_tested": grp["species"].nunique(),
                        "n_species_fdr05": int(grp["fdr"].lt(0.05).sum()),
                        "species": ",".join(sorted(grp["species"])),
                        "mean_fold": float(grp["fold_enrichment"].mean()),
                        "min_fdr": float(grp["fdr"].min()),
                        "combined_p": partial_conjunction(grp["p_value"].tolist()),
                    }
                )

            # theme-level enrichment
            fg_genes = {
                species: set(
                    members[
                        members["species"].eq(species)
                        & members["direction"].eq(direction)
                        & members["evidence_mode"].eq(mode)
                    ]["gene_id"]
                    .str.split("__", n=1)
                    .str[-1]
                )
                for species in SPECIES
            }
            for species in SPECIES:
                go = go_by_species[species]
                all_genes = set(go["protein_id"])
                fg = fg_genes[species] & all_genes
                if not fg:
                    continue
                for theme, keywords in THEMES.items():
                    theme_go = go[
                        go["go_description"]
                        .str.lower()
                        .apply(lambda d: any(k in d for k in keywords))
                    ]
                    theme_genes = set(theme_go["protein_id"]) & all_genes
                    if not theme_genes:
                        continue
                    n, N = len(fg), len(all_genes)
                    k = len(fg & theme_genes)
                    K = len(theme_genes)
                    theme_rows.append(
                        {
                            "evidence_mode": mode,
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
        meta["fdr"] = meta.groupby(["evidence_mode", "direction"])[
            "combined_p"
        ].transform(bh_group)
        meta.sort_values(["evidence_mode", "direction", "combined_p"]).to_csv(
            OUT / "GO_nonoverlap_meta.tsv", sep="\t", index=False
        )

    themes = pd.DataFrame(theme_rows)
    if not themes.empty:
        # Each species tests the seven themes within each direction/evidence mode.
        themes["fdr"] = themes.groupby(["evidence_mode", "direction", "species"])[
            "p_value"
        ].transform(bh_group)
        meta_theme = (
            themes.groupby(["evidence_mode", "direction", "theme"])
            .agg(
                n_species_tested=("species", "nunique"),
                n_species_fdr05=("fdr", lambda v: int((v < 0.05).sum())),
                mean_fold=("fold_enrichment", "mean"),
                min_fdr=("fdr", "min"),
                combined_p=("p_value", lambda v: partial_conjunction(list(v))),
            )
            .reset_index()
        )
        meta_theme["fdr_combined"] = meta_theme.groupby(["evidence_mode", "direction"])[
            "combined_p"
        ].transform(bh_group)
        themes.to_csv(OUT / "GO_nonoverlap_theme.tsv", sep="\t", index=False)
        meta_theme.to_csv(OUT / "GO_nonoverlap_theme_meta.tsv", sep="\t", index=False)

    if not meta.empty:
        print("== term-level meta (all evidence) ==")
        m = meta[meta["evidence_mode"].eq("all")].sort_values("combined_p")
        print(
            m.head(15)[
                [
                    "direction",
                    "go_description",
                    "n_species_tested",
                    "n_species_fdr05",
                    "mean_fold",
                    "combined_p",
                    "fdr",
                ]
            ].to_string(index=False)
        )
    if "meta_theme" in dir():
        print()
        print("== theme-level (all evidence) ==")
        t = meta_theme[meta_theme["evidence_mode"].eq("all")].sort_values("combined_p")
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
