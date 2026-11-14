#!/usr/bin/env python3
"""KEGG pathway enrichment of genes in non-overlapping A/B retention regions.

Mirrors the GO workflow (regional_go.py): regions are defined on the
FraBias-style curves (every overlapping 100-gene window along the GTR target,
collinear anchors). A window is non-overlapping when |A - B| > 5 percentage
points; genes are assigned the class of the window centred on them. For
A-favoured regions the A-copy genes are tested, for B-favoured regions the
B-copy genes. Background: all pathway-annotated genes of the species (whole
genome), plus an anchored-subgenome background as a second scope.

Annotation: KEGG gene IDs transferred from Swiss-Prot cross-references,
resolved to KO / reference pathways by kegg_annotations.py.
"""

from __future__ import annotations

from enrichment_utils import bh_group, fisher_greater, hypergeom_sf

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
ROOT = args.root
OUT = ROOT / "GO"
BASE = (
    ROOT
    / "downstream_analysis"
    / "figure2_polyploid"
    / "strict_wgd_mapping"
    / "04_fractbias_8track_reference_exact"
)
ORDER = BASE / "01_catalog" / "gtr_gene_order.tsv"
NOLOCAL = BASE / "02_last" / "GTR.nolocaldups.bed"
QBED = BASE / "02_last" / "query8.nolocaldups.bed"
MERGED = BASE / "04_quota_align" / "query8_vs_GTR.initial.qa.merged"
ANN = ROOT / "annotation" / "functional"

SPECIES = ["PZI", "PSP", "PPA"]
SHORT = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
CHROMS = [f"GTR__Chr{i:02d}" for i in range(1, 10)]
WINDOW = 100
HALF = WINDOW // 2
DIVERGENT = 0.05


def read_pathways(species: str) -> pd.DataFrame:
    table = pd.read_csv(
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
    )
    return table.dropna(subset=["protein_id", "pathway_id"]).drop_duplicates(
        ["protein_id", "pathway_id"]
    )


def main() -> None:
    order = pd.read_csv(ORDER, sep="\t")
    gene2gi = dict(zip(order["gtr_gene"], order["gtr_global_index"]))
    bed = pd.read_csv(
        NOLOCAL,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    tid2gi = bed["gene"].map(gene2gi).to_numpy()
    qbed = pd.read_csv(
        QBED,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "gene"],
        usecols=[0, 1, 2, 3],
    )
    qgene = qbed["gene"].to_numpy()

    merged = pd.read_csv(
        MERGED,
        sep="\t",
        skiprows=1,
        header=None,
        usecols=[1, 2, 3],
        names=["tidx", "query", "qidx"],
    ).dropna(subset=["query"])
    merged["track"] = merged["query"].str.split("__").str[0]
    merged = merged[merged["track"].str.match(r"^(PZI|PSP|PPA)_[AB]$")].copy()
    merged["gi"] = (
        merged["tidx"]
        .astype(int)
        .map(
            lambda t: (
                int(tid2gi[t])
                if 0 <= t < len(tid2gi) and not pd.isna(tid2gi[t])
                else np.nan
            )
        )
    )
    merged["gene_id"] = (
        merged["qidx"]
        .astype(int)
        .map(lambda i: qgene[i].split("__", 1)[-1] if 0 <= i < len(qgene) else np.nan)
    )
    merged = merged.dropna(subset=["gi", "gene_id"])
    merged["gi"] = merged["gi"].astype(int)

    presence = {}
    for species in SPECIES:
        for copy in ("A", "B"):
            track = f"{species}_{copy}"
            presence[track] = set(merged.loc[merged["track"].eq(track), "gi"])

    class_of_gene: dict[str, np.ndarray] = {}
    for species in SPECIES:
        cls = np.full(len(order), "none", dtype=object)
        for chrom in CHROMS:
            sub = order[order["gtr_seqid"].eq(chrom)].sort_values("gtr_global_index")
            idx = sub.index.to_numpy()
            gi_list = sub["gtr_global_index"].to_numpy()
            n = len(gi_list)
            if n < WINDOW:
                continue
            a = np.isin(gi_list, list(presence[f"{species}_A"])).astype(float)
            b = np.isin(gi_list, list(presence[f"{species}_B"])).astype(float)
            ca = np.concatenate([[0], np.cumsum(a)])
            cb = np.concatenate([[0], np.cumsum(b)])
            wa = (ca[WINDOW:] - ca[:-WINDOW]) / WINDOW
            wb = (cb[WINDOW:] - cb[:-WINDOW]) / WINDOW
            nwin = n - WINDOW + 1
            for i in range(n):
                s = i - HALF
                if s < 0 or s >= nwin:
                    continue
                diff = wa[s] - wb[s]
                if diff > DIVERGENT:
                    cls[idx[i]] = "A_favoured"
                elif -diff > DIVERGENT:
                    cls[idx[i]] = "B_favoured"
        class_of_gene[species] = cls

    rows_ora: list = []
    rows_contrast: list = []
    rows_summary: list = []
    rows_members: list = []
    pathway_by_species = {}
    meta = {}
    for species in SPECIES:
        pathways = read_pathways(SHORT[species])
        pathway_by_species[species] = (
            pathways.groupby("protein_id")["pathway_id"].agg(set).to_dict()
        )
        for r in pathways.drop_duplicates("pathway_id").itertuples(index=False):
            meta[r.pathway_id] = (r.pathway_name, r.pathway_class_a, r.pathway_class_b)

    for species in SPECIES:
        sub = merged[merged["track"].str.startswith(species)].copy()
        sub["cls"] = class_of_gene[species][sub["gi"].to_numpy()]
        for direction, fav_copy, other_copy in (
            ("A_favoured", "A", "B"),
            ("B_favoured", "B", "A"),
        ):
            fav = sub[
                sub["track"].eq(f"{species}_{fav_copy}") & sub["cls"].eq(direction)
            ]
            other = sub[
                sub["track"].eq(f"{species}_{other_copy}") & sub["cls"].eq(direction)
            ]
            shared_ids = set(fav["gi"]) & set(other["gi"])
            fg_all = set(fav["gene_id"]) & set(pathway_by_species[species])
            fg_shared = set(fav.loc[fav["gi"].isin(shared_ids), "gene_id"]) & set(
                pathway_by_species[species]
            )
            bg_genome = set(pathway_by_species[species])
            bg_anchored = {
                qgene[i].split("__", 1)[-1]
                for i in merged.loc[
                    merged["track"].str.startswith(species), "qidx"
                ].astype(int)
                if 0 <= i < len(qgene)
            } & set(pathway_by_species[species])
            meta_sets = {"all_anchored": fg_all, "shared_only": fg_shared}
            for scope, fg in meta_sets.items():
                for bg_name, bg in (("genome", bg_genome), ("anchored", bg_anchored)):
                    units = {g: pathway_by_species[species].get(g, set()) for g in bg}
                    index: dict[str, set] = defaultdict(set)
                    for g, terms in units.items():
                        for t in terms:
                            index[t].add(g)
                    n = len(fg)
                    N = len(bg)
                    if n == 0:
                        continue
                    for pathway_id, ids in index.items():
                        k = len(fg & ids)
                        K = len(ids)
                        name, class_a, class_b = meta.get(pathway_id, ("", "", ""))
                        rows_ora.append(
                            {
                                "species": species,
                                "direction": direction,
                                "scope": scope,
                                "background": bg_name,
                                "pathway_id": pathway_id,
                                "pathway_name": name,
                                "pathway_class_a": class_a,
                                "pathway_class_b": class_b,
                                "foreground_genes": n,
                                "background_genes": N,
                                "term_foreground": k,
                                "term_background": K,
                                "fold_enrichment": (k / n) / (K / N),
                                "p_value": hypergeom_sf(k, N, K, n),
                            }
                        )
            rows_summary.append(
                {
                    "species": species,
                    "direction": direction,
                    "foreground_genes_with_pathway": len(fg_all),
                    "foreground_shared_with_pathway": len(fg_shared),
                    "genome_background_with_pathway": len(bg_genome),
                    "anchored_background_with_pathway": len(bg_anchored),
                }
            )
            members = fav[["track", "gi", "gene_id", "cls"]].drop_duplicates()
            members["species"] = species
            members["direction"] = direction
            rows_members.append(members)

        for fav_copy, direction in (("A", "A_favoured"), ("B", "B_favoured")):
            fav = sub[
                sub["track"].eq(f"{species}_{fav_copy}") & sub["cls"].eq(direction)
            ]
            genes = set(fav["gene_id"])
            if fav_copy == "A":
                a_genes = genes
            else:
                b_genes = genes
        units = {
            g: pathway_by_species[species].get(g, set()) for g in (a_genes | b_genes)
        }
        index: dict[str, set] = defaultdict(set)
        for g, terms in units.items():
            for t in terms:
                index[t].add(g)
        for pathway_id, ids in index.items():
            a = len(a_genes & ids)
            c = len(b_genes & ids)
            if a + c == 0:
                continue
            b = len(a_genes) - a
            d = len(b_genes) - c
            name, class_a, class_b = meta.get(pathway_id, ("", "", ""))
            rows_contrast.append(
                {
                    "species": species,
                    "pathway_id": pathway_id,
                    "pathway_name": name,
                    "pathway_class_a": class_a,
                    "pathway_class_b": class_b,
                    "A_favoured_genes": len(a_genes),
                    "B_favoured_genes": len(b_genes),
                    "term_A_favoured": a,
                    "term_B_favoured": c,
                    "odds_ratio_A_over_B": ((a + 0.5) * (d + 0.5))
                    / ((b + 0.5) * (c + 0.5)),
                    "p_A_enriched": fisher_greater(a, b, c, d),
                    "p_B_enriched": fisher_greater(c, d, a, b),
                    "p_two_sided": min(
                        1.0,
                        2.0
                        * min(fisher_greater(a, b, c, d), fisher_greater(c, d, a, b)),
                    ),
                }
            )

    ora = pd.DataFrame(rows_ora)
    if not ora.empty:
        ora["fdr"] = ora.groupby(
            ["species", "direction", "scope", "background"],
            group_keys=False,
        )["p_value"].transform(bh_group)
        ora.sort_values(["species", "direction", "background", "fdr"]).to_csv(
            OUT / "KEGG_nonoverlap_ORA.tsv", sep="\t", index=False
        )
    contrast = pd.DataFrame(rows_contrast)
    if not contrast.empty:
        contrast["log2_odds_ratio_A_over_B"] = np.log2(contrast["odds_ratio_A_over_B"])
        contrast["fdr_two_sided"] = contrast.groupby(["species"], group_keys=False)[
            "p_two_sided"
        ].transform(bh_group)
        contrast.sort_values(["species", "fdr_two_sided"]).to_csv(
            OUT / "KEGG_nonoverlap_contrast.tsv", sep="\t", index=False
        )
    pd.DataFrame(rows_summary).to_csv(
        OUT / "KEGG_nonoverlap_summary.tsv", sep="\t", index=False
    )
    if rows_members:
        pd.concat(rows_members, ignore_index=True).to_csv(
            OUT / "KEGG_nonoverlap_gene_members.tsv", sep="\t", index=False
        )

    print(pd.DataFrame(rows_summary).to_string(index=False))
    if not ora.empty:
        print()
        print("ORA tests:", len(ora), "| FDR<0.05:", int(ora.fdr.lt(0.05).sum()))
        sig = ora[ora.fdr.lt(0.05)].sort_values("fdr")
        if not sig.empty:
            print(
                sig.head(25)[
                    [
                        "species",
                        "direction",
                        "scope",
                        "background",
                        "pathway_name",
                        "term_foreground",
                        "foreground_genes",
                        "term_background",
                        "background_genes",
                        "fold_enrichment",
                        "fdr",
                    ]
                ].to_string(index=False)
            )
    if not contrast.empty:
        print()
        print(
            "contrast tests:",
            len(contrast),
            "| FDR<0.05:",
            int(contrast.fdr_two_sided.lt(0.05).sum()),
        )
        print(
            contrast.sort_values("p_two_sided")
            .head(10)[
                [
                    "species",
                    "pathway_name",
                    "term_A_favoured",
                    "term_B_favoured",
                    "odds_ratio_A_over_B",
                    "p_two_sided",
                    "fdr_two_sided",
                ]
            ]
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
