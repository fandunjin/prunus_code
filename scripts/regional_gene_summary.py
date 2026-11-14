#!/usr/bin/env python3
"""Extract genes in non-overlapping (|A-B| > 5 pp) regions, split by direction.

Directions follow regional_go.py: A_favoured = A > B, B_favoured = B > A.
Gene assignment replicates the regional_go.py window logic exactly
(100-gene overlapping windows on the collinear-anchor curves; each target gene
inherits the class of the window centred on its position, offset HALF = 50).

Outputs (all under GO/0915/):
  genes_A_gt_B.tsv / genes_B_gt_A.tsv
      all anchored query genes of the three species in the divergent windows,
      both copies, flagged favoured_track, with per-window retention values,
      delta (A - B, pp) and swissprot annotation columns.
  by_species/{PZI,PSP,PPA}_{A_gt_B,B_gt_A}.tsv
      favoured-track (foreground-equivalent) gene sets, one file per species.
  annotation_summary.tsv
      counts and annotation coverage per species x direction.
  GO_ORA_*.tsv / KEGG_ORA_*.tsv
      subsets of the regenerated ORA tables for each direction.
  top_terms_compact.tsv / top_term_genes.tsv
      top enriched GO/KEGG terms and their member genes.
  SUMMARY.md
      short human-readable digest.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
ROOT = args.root
GO = ROOT / "GO"
OUT = GO / "0915"
ANN = ROOT / "annotation" / "functional"
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

SPECIES = ["PZI", "PSP", "PPA"]
SHORT = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
CHROMS = [f"GTR__Chr{i:02d}" for i in range(1, 10)]
WINDOW = 100
HALF = WINDOW // 2
DIVERGENT = 0.05

DIRECTION_LABEL = {"A_favoured": "A_gt_B", "B_favoured": "B_gt_A"}


def load_annotation(species: str):
    short = SHORT[species]
    func = pd.read_csv(
        ANN / short / "swissprot_functional_annotation.tsv",
        sep="\t",
        usecols=["protein_id", "swissprot_entry", "gene_names", "protein_description"],
        dtype=str,
    ).drop_duplicates("protein_id")
    go = pd.read_csv(
        ANN / short / "swissprot_go.tsv",
        sep="\t",
        usecols=["protein_id", "go_id"],
        dtype=str,
    ).drop_duplicates(["protein_id", "go_id"])
    go_n = go.groupby("protein_id")["go_id"].nunique().rename("go_terms_n")
    kegg = pd.read_csv(
        ANN / short / "swissprot_kegg_pathway.tsv",
        sep="\t",
        usecols=["protein_id", "pathway_id", "pathway_name"],
        dtype=str,
    ).drop_duplicates(["protein_id", "pathway_id"])
    kegg_n = (
        kegg.groupby("protein_id")["pathway_id"].nunique().rename("kegg_pathways_n")
    )
    kegg_top = (
        kegg.drop_duplicates("protein_id")
        .set_index("protein_id")["pathway_name"]
        .rename("kegg_top_pathway")
    )
    return func, go, go_n, kegg, kegg_n, kegg_top


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "by_species").mkdir(exist_ok=True)

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
    merged["copy"] = merged["track"].str[-1]

    presence = {}
    for species in SPECIES:
        for copy in ("A", "B"):
            track = f"{species}_{copy}"
            presence[track] = set(merged.loc[merged["track"].eq(track), "gi"])

    cls_of_gene: dict[str, np.ndarray] = {}
    delta_of_gene: dict[str, np.ndarray] = {}
    wa_of_gene: dict[str, np.ndarray] = {}
    wb_of_gene: dict[str, np.ndarray] = {}
    for species in SPECIES:
        cls = np.full(len(order), "none", dtype=object)
        delta = np.full(len(order), np.nan)
        wa_g = np.full(len(order), np.nan)
        wb_g = np.full(len(order), np.nan)
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
                else:
                    continue
                delta[idx[i]] = diff * 100.0
                wa_g[idx[i]] = wa[s] * 100.0
                wb_g[idx[i]] = wb[s] * 100.0
        cls_of_gene[species] = cls
        delta_of_gene[species] = delta
        wa_of_gene[species] = wa_g
        wb_of_gene[species] = wb_g

    rows = []
    for species in SPECIES:
        sub = merged[merged["track"].str.startswith(species)].copy()
        sub["window_class"] = cls_of_gene[species][sub["gi"].to_numpy()]
        sub["window_delta_pp"] = delta_of_gene[species][sub["gi"].to_numpy()]
        sub["window_A_pct"] = wa_of_gene[species][sub["gi"].to_numpy()]
        sub["window_B_pct"] = wb_of_gene[species][sub["gi"].to_numpy()]
        sub["species"] = species
        sub = sub[sub["window_class"].ne("none")]
        rows.append(sub)
    genes = pd.concat(rows, ignore_index=True)
    genes = genes.rename(columns={"gi": "gtr_gi"})

    frames = {}
    for direction in ("A_favoured", "B_favoured"):
        d = genes[genes["window_class"].eq(direction)].copy()
        d["direction"] = DIRECTION_LABEL[direction]
        d["favoured_track"] = d["copy"].eq("A" if direction == "A_favoured" else "B")
        d["abs_delta"] = d["window_delta_pp"].abs()
        d = (
            d.sort_values("abs_delta", ascending=False)
            .drop_duplicates(subset=["species", "track", "gene_id"])
            .drop(columns=["abs_delta"])
        )
        frames[direction] = d

    func_all, go_all, go_n_all, kegg_all, kegg_n_all, kegg_top_all = (
        {},
        {},
        {},
        {},
        {},
        {},
    )
    for species in SPECIES:
        func, go, go_n, kegg, kegg_n, kegg_top = load_annotation(species)
        func_all[species], go_all[species] = func, go
        go_n_all[species], kegg_all[species] = go_n, kegg
        kegg_n_all[species], kegg_top_all[species] = kegg_n, kegg_top

    def annotate(df: pd.DataFrame) -> pd.DataFrame:
        out = []
        for species, sub in df.groupby("species"):
            sub = sub.merge(
                func_all[species], left_on="gene_id", right_on="protein_id", how="left"
            )
            sub = sub.merge(
                go_n_all[species], left_on="gene_id", right_index=True, how="left"
            )
            sub = sub.merge(
                kegg_n_all[species], left_on="gene_id", right_index=True, how="left"
            )
            sub = sub.merge(
                kegg_top_all[species], left_on="gene_id", right_index=True, how="left"
            )
            sub["go_terms_n"] = sub["go_terms_n"].fillna(0).astype(int)
            sub["kegg_pathways_n"] = sub["kegg_pathways_n"].fillna(0).astype(int)
            out.append(sub)
        return pd.concat(out, ignore_index=True)

    cols = [
        "species",
        "direction",
        "track",
        "copy",
        "favoured_track",
        "gene_id",
        "gtr_gi",
        "window_delta_pp",
        "window_A_pct",
        "window_B_pct",
        "swissprot_entry",
        "gene_names",
        "protein_description",
        "go_terms_n",
        "kegg_pathways_n",
        "kegg_top_pathway",
    ]

    annotated = {}
    for direction, label in (("A_favoured", "A_gt_B"), ("B_favoured", "B_gt_A")):
        d = annotate(frames[direction])[cols].copy()
        d["protein_description"] = d["protein_description"].str.slice(0, 160)
        sort_asc = label == "B_gt_A"
        d = d.sort_values(
            ["species", "favoured_track", "window_delta_pp"],
            ascending=[True, False, sort_asc],
        )
        d.to_csv(OUT / f"genes_{label}.tsv", sep="\t", index=False)
        annotated[direction] = d
        fav = d[d["favoured_track"]].copy()
        fav.to_csv(OUT / f"genes_foreground_{label}.tsv", sep="\t", index=False)
        for species in SPECIES:
            fav[fav["species"].eq(species)].to_csv(
                OUT / "by_species" / f"{species}_{label}.tsv", sep="\t", index=False
            )

    summary_rows = []
    for direction, label in (("A_favoured", "A_gt_B"), ("B_favoured", "B_gt_A")):
        d = annotated[direction]
        for species in SPECIES:
            sub = d[d["species"].eq(species)]
            fav = sub[sub["favoured_track"]]
            summary_rows.append(
                {
                    "species": species,
                    "direction": label,
                    "n_genes_all_copies": len(sub),
                    "n_genes_favoured_track": len(fav),
                    "n_other_copy_same_windows": len(sub) - len(fav),
                    "n_favoured_annotated": int(fav["swissprot_entry"].notna().sum()),
                    "pct_favoured_annotated": (
                        round(100.0 * fav["swissprot_entry"].notna().mean(), 1)
                        if len(fav)
                        else np.nan
                    ),
                    "n_favoured_with_GO": int((fav["go_terms_n"] > 0).sum()),
                    "n_favoured_with_KEGG": int((fav["kegg_pathways_n"] > 0).sum()),
                    "mean_abs_delta_pp_favoured": (
                        round(float(fav["window_delta_pp"].abs().mean()), 2)
                        if len(fav)
                        else np.nan
                    ),
                    "max_abs_delta_pp_favoured": (
                        round(float(fav["window_delta_pp"].abs().max()), 2)
                        if len(fav)
                        else np.nan
                    ),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "annotation_summary.tsv", sep="\t", index=False)

    go_ora = pd.read_csv(GO / "GO_nonoverlap_ORA.tsv", sep="\t")
    kegg_ora = pd.read_csv(GO / "KEGG_nonoverlap_ORA.tsv", sep="\t")
    go_ora[go_ora["direction"].eq("A_favoured")].sort_values(
        ["evidence_mode", "species", "background", "fdr"]
    ).to_csv(OUT / "GO_ORA_A_gt_B.tsv", sep="\t", index=False)
    go_ora[go_ora["direction"].eq("B_favoured")].sort_values(
        ["evidence_mode", "species", "background", "fdr"]
    ).to_csv(OUT / "GO_ORA_B_gt_A.tsv", sep="\t", index=False)
    kegg_ora[kegg_ora["direction"].eq("A_favoured")].sort_values(
        ["species", "background", "fdr"]
    ).to_csv(OUT / "KEGG_ORA_A_gt_B.tsv", sep="\t", index=False)
    kegg_ora[kegg_ora["direction"].eq("B_favoured")].sort_values(
        ["species", "background", "fdr"]
    ).to_csv(OUT / "KEGG_ORA_B_gt_A.tsv", sep="\t", index=False)

    compact = []
    go_sel = go_ora[
        go_ora["evidence_mode"].eq("all")
        & go_ora["scope"].eq("all_anchored")
        & go_ora["background"].eq("anchored")
    ]
    for (direction, species), sub in go_sel.groupby(["direction", "species"]):
        for rank, (_, r) in enumerate(sub.sort_values("fdr").head(8).iterrows(), 1):
            compact.append(
                {
                    "domain": "GO",
                    "direction": DIRECTION_LABEL[direction],
                    "species": species,
                    "rank": rank,
                    "term_id": r["go_id"],
                    "term": r["go_description"],
                    "aspect_class": r["go_aspect"],
                    "term_genes": r["term_foreground"],
                    "fold_enrichment": round(r["fold_enrichment"], 2),
                    "fdr": f"{r['fdr']:.3g}",
                }
            )
    kegg_sel = kegg_ora[
        kegg_ora["background"].eq("anchored") & kegg_ora["scope"].eq("all_anchored")
    ]
    for (direction, species), sub in kegg_sel.groupby(["direction", "species"]):
        for rank, (_, r) in enumerate(sub.sort_values("fdr").head(5).iterrows(), 1):
            compact.append(
                {
                    "domain": "KEGG",
                    "direction": DIRECTION_LABEL[direction],
                    "species": species,
                    "rank": rank,
                    "term_id": r["pathway_id"],
                    "term": r["pathway_name"],
                    "aspect_class": r["pathway_class_b"],
                    "term_genes": r["term_foreground"],
                    "fold_enrichment": round(r["fold_enrichment"], 2),
                    "fdr": f"{r['fdr']:.3g}",
                }
            )
    compact = pd.DataFrame(compact)
    compact.to_csv(OUT / "top_terms_compact.tsv", sep="\t", index=False)

    gene_rows = []
    for direction, label in (("A_favoured", "A_gt_B"), ("B_favoured", "B_gt_A")):
        d = annotated[direction]
        fav = d[d["favoured_track"]]
        for species in SPECIES:
            terms = compact[
                compact["domain"].eq("GO")
                & compact["direction"].eq(label)
                & compact["species"].eq(species)
                & compact["rank"].le(5)
            ]
            go_map = go_all[species].groupby("go_id")["protein_id"].agg(set).to_dict()
            for _, t in terms.iterrows():
                members = go_map.get(t["term_id"], set()) & set(
                    fav.loc[fav["species"].eq(species), "gene_id"]
                )
                sub = fav[fav["gene_id"].isin(members)]
                sub = sub.sort_values(
                    "window_delta_pp", ascending=(label == "B_gt_A")
                ).head(12)
                for _, g in sub.iterrows():
                    gene_rows.append(
                        {
                            "direction": label,
                            "species": species,
                            "GO_term": t["term"],
                            "GO_id": t["term_id"],
                            "gene_id": g["gene_id"],
                            "gene_description": (
                                g["protein_description"]
                                if isinstance(g["protein_description"], str)
                                else ""
                            )[:90],
                            "window_delta_pp": round(g["window_delta_pp"], 1),
                        }
                    )
    pd.DataFrame(gene_rows).to_csv(OUT / "top_term_genes.tsv", sep="\t", index=False)

    theme_go = pd.read_csv(GO / "GO_nonoverlap_theme_meta.tsv", sep="\t")
    theme_go["direction_label"] = theme_go["direction"].map(DIRECTION_LABEL)
    theme_go.sort_values(["direction", "fdr_combined"]).to_csv(
        OUT / "theme_meta_GO.tsv", sep="\t", index=False
    )
    theme_kegg = pd.read_csv(GO / "KEGG_nonoverlap_theme_meta.tsv", sep="\t")
    theme_kegg["direction_label"] = theme_kegg["direction"].map(DIRECTION_LABEL)
    theme_kegg.sort_values(["direction", "fdr_combined"]).to_csv(
        OUT / "theme_meta_KEGG.tsv", sep="\t", index=False
    )

    lines = ["# Non-overlapping (|A-B| > 5 pp) region genes, split by direction", ""]
    lines.append(
        "A_gt_B = A_favoured (A track higher), B_gt_A = B_favoured (B track higher)."
    )
    lines.append(
        "Threshold: |A-B| > 5 percentage points, 100-gene collinear-anchor windows,"
    )
    lines.append(
        "gene class = window centred on the gene's GTR position (same logic as regional_go.py)."
    )
    lines.append("")
    lines.append("## Counts and annotation coverage (favoured track)")
    lines.append("")
    lines.append(summary.to_markdown(index=False))
    lines.append("")
    lines.append("## Top terms (scope all_anchored, background anchored)")
    lines.append("")
    lines.append(compact.head(40).to_markdown(index=False))
    lines.append("")
    theme_all = theme_go[theme_go["evidence_mode"].eq("all")].sort_values(
        ["direction", "fdr_combined"]
    )
    lines.append("## Cross-species theme meta (Fisher, species evidence_mode = all)")
    lines.append("")
    lines.append(
        theme_all[
            ["direction_label", "theme", "n_species_fdr05", "mean_fold", "fdr_combined"]
        ].to_markdown(index=False)
    )
    lines.append("")
    lines.append("## Files")
    lines.append(
        "- `genes_A_gt_B.tsv` / `genes_B_gt_A.tsv`: all anchored genes (both copies) in divergent windows + annotation"
    )
    lines.append(
        "- `genes_foreground_*.tsv`: favoured-track genes only (regional_go.py foreground equivalent)"
    )
    lines.append("- `by_species/*.tsv`: per-species favoured-track gene lists")
    lines.append("- `annotation_summary.tsv`: counts per species x direction")
    lines.append(
        "- `GO_ORA_*.tsv` / `KEGG_ORA_*.tsv`: full ORA subsets from the regenerated pipeline"
    )
    lines.append("- `top_terms_compact.tsv`: top GO/KEGG terms per species x direction")
    lines.append("- `top_term_genes.tsv`: member genes of the top-5 GO terms per group")
    lines.append(
        "- `theme_meta_GO.tsv` / `theme_meta_KEGG.tsv`: cross-species theme-level meta (same source as Fig 3e)"
    )
    (OUT / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")

    print(summary.to_string(index=False))
    print()
    print("files written to", OUT)


if __name__ == "__main__":
    main()
