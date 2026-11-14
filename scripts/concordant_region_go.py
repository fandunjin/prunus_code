#!/usr/bin/env python3
"""GO enrichment for Chr03+Chr07 regions with concordant A/B retention trends.

Selection: for every 10-gene bin on GTR Chr03/Chr07, the local slope of the
11-bin-smoothed A and B retention curves is estimated over a +/-5-bin
neighbourhood (least squares). Bins where both slopes are positive (concordant
up / co-retention) or both negative (concordant down / co-loss) are selected,
per species and without requiring cross-species consensus. Genes of those bins
are tested per species (one-sided hypergeometric, BH) against two backgrounds
(whole-genome GO genes and Chr03+Chr07 anchored genes), for scopes
shared A/B anchors and any anchored copy, and combined across species with
a Bonferroni partial-conjunction test (at least two of three species).

Source data: GO/window_retention_profile.tsv, figures/figure4_panel_a_source_data.tsv,
downstream_analysis/.../quota_filtered_anchors.tsv, annotation/functional/*/swissprot_go.tsv
"""

from __future__ import annotations

from enrichment_utils import bh_group, hypergeom_sf, partial_conjunction

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
ANCHORS = BASE / "04_quota_align" / "normalized" / "quota_filtered_anchors.tsv"
WINDOWS = OUT / "window_retention_profile.tsv"
SOURCE = ROOT / "figures" / "figure4_panel_a_source_data.tsv"
GO_DIR = ROOT / "annotation" / "functional"

SPECIES = ["PZI", "PSP", "PPA"]
SHORT = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
EVIDENCE_MODES = ("all", "non_IEA")
CHROMS = ["GTR__Chr03", "GTR__Chr07"]
TREND_RADIUS = 5
MIN_NEIGHBOUR_BINS = 6
MIN_SLOPE_CHANGE = 0.01
MIN_TERM_BG = 10
MIN_SPECIES = 2
STRATA = ["concordant_up", "concordant_down", "concordant_any"]
SCOPES = ["shared", "any"]


def read_go(path: Path, mode: str) -> pd.DataFrame:
    usecols = ["protein_id", "go_id", "go_aspect", "go_description", "go_evidence"]
    go = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str)
    go = go.dropna(subset=["protein_id", "go_id", "go_aspect"])
    if mode == "non_IEA":
        go = go[go["go_evidence"].fillna("").str.upper().ne("IEA")]
    return go.drop_duplicates(["protein_id", "go_id"])


def slope(values: np.ndarray) -> float:
    mask = ~np.isnan(values)
    if mask.sum() < MIN_NEIGHBOUR_BINS:
        return float("nan")
    x = np.arange(len(values), dtype=float)[mask]
    y = values[mask]
    x = x - x.mean()
    return float((x * (y - y.mean())).sum() / (x * x).sum())


def trend_table(windows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (chrom, species), sub in windows.groupby(["chromosome", "species"]):
        sub = sub.sort_values("bin").reset_index(drop=True)
        a = sub["A_smoothed"].to_numpy(dtype=float)
        b = sub["B_smoothed"].to_numpy(dtype=float)
        elig = sub["eligible"].to_numpy(dtype=bool)
        n = len(sub)
        slope_a = np.full(n, np.nan)
        slope_b = np.full(n, np.nan)
        for i in range(n):
            lo, hi = i - TREND_RADIUS, i + TREND_RADIUS + 1
            lo_c, hi_c = max(0, lo), min(n, hi)
            aa = np.where(elig[lo_c:hi_c], a[lo_c:hi_c], np.nan)
            bb = np.where(elig[lo_c:hi_c], b[lo_c:hi_c], np.nan)
            slope_a[i] = slope(aa)
            slope_b[i] = slope(bb)
        for i in range(n):
            if not elig[i] or np.isnan(slope_a[i]) or np.isnan(slope_b[i]):
                cls = "excluded"
            elif (
                slope_a[i] * TREND_RADIUS * 2 >= MIN_SLOPE_CHANGE
                and slope_b[i] * TREND_RADIUS * 2 >= MIN_SLOPE_CHANGE
            ):
                cls = "concordant_up"
            elif (
                slope_a[i] * TREND_RADIUS * 2 <= -MIN_SLOPE_CHANGE
                and slope_b[i] * TREND_RADIUS * 2 <= -MIN_SLOPE_CHANGE
            ):
                cls = "concordant_down"
            else:
                cls = "inconsistent"
            rows.append(
                {
                    "chromosome": chrom,
                    "species": species,
                    "bin": int(sub.at[i, "bin"]),
                    "pillar_start": int(sub.at[i, "pillar_start"]),
                    "pillar_end": int(sub.at[i, "pillar_end"]),
                    "A_smoothed": float(sub.at[i, "A_smoothed"]),
                    "B_smoothed": float(sub.at[i, "B_smoothed"]),
                    "slope_A": slope_a[i],
                    "slope_B": slope_b[i],
                    "trend_class": cls,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    windows = pd.read_csv(WINDOWS, sep="\t")
    windows = windows[windows["chromosome"].isin(CHROMS)].copy()
    trends = trend_table(windows)
    trends.to_csv(OUT / "GO_chr3_chr7_concordant_windows.tsv", sep="\t", index=False)

    source = pd.read_csv(SOURCE, sep="\t")

    anchors = pd.read_csv(ANCHORS, sep="\t", dtype=str)
    anchors["gtr_global_index"] = pd.to_numeric(
        anchors["gtr_global_index"], errors="coerce"
    )
    anchors = anchors[
        anchors["gtr_seqid"].isin(CHROMS)
        & anchors["track_id"].str.match(r"^(PZI|PSP|PPA)_[AB]$")
    ].dropna(subset=["gtr_global_index", "query_gene"])
    anchors["gtr_global_index"] = anchors["gtr_global_index"].astype(int)
    anchors["species"] = anchors["track_id"].str.split("_").str[0]
    anchors["copy"] = anchors["track_id"].str.split("_").str[1]
    anchors["gene_id"] = anchors["query_gene"].str.split("__", n=1).str[-1]

    copy_index: dict[str, dict[int, dict[str, set[str]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(set))
    )
    for r in anchors.itertuples(index=False):
        copy_index[r.species][r.gtr_global_index][r.copy].add(r.gene_id)

    bin_pillars: dict[str, list[int]] = {}
    for r in trends.itertuples(index=False):
        ids = source.iloc[r.pillar_start : r.pillar_end]["gtr_global_index"].astype(int)
        bin_pillars[(r.species, r.bin)] = list(ids)

    strata_members: dict[str, dict[tuple[str, str], set[int]]] = {}
    for species in SPECIES:
        bins = trends[trends["species"].eq(species)]
        up_bins = set(bins.loc[bins["trend_class"].eq("concordant_up"), "bin"])
        down_bins = set(bins.loc[bins["trend_class"].eq("concordant_down"), "bin"])
        up = (
            set().union(*[bin_pillars[(species, b)] for b in up_bins])
            if up_bins
            else set()
        )
        down = (
            set().union(*[bin_pillars[(species, b)] for b in down_bins])
            if down_bins
            else set()
        )
        strata_members[species] = {
            "concordant_up": up,
            "concordant_down": down,
            "concordant_any": up | down,
        }

    ora_rows: list = []
    member_rows: list = []
    summary_rows: list = []
    for mode in EVIDENCE_MODES:
        go_by_species = {}
        meta: dict[str, tuple[str, str]] = {}
        for species in SPECIES:
            go = read_go(GO_DIR / SHORT[species] / "swissprot_go.tsv", mode)
            go_by_species[species] = (
                go.groupby("protein_id")["go_id"].agg(set).to_dict()
            )
            for r in go.itertuples(index=False):
                meta[r.go_id] = (r.go_aspect, r.go_description)

        per_species_results: dict[tuple[str, str, str], dict[str, set[str]]] = {}
        for species in SPECIES:
            genome_bg = set(go_by_species[species])
            shared_pillars = {
                gi
                for gi, copies in copy_index[species].items()
                if "A" in copies and "B" in copies
            }
            any_pillars = set(copy_index[species])
            backgrounds = {
                "genome": genome_bg,
                "chr37_anchored": None,
            }
            for stratum in STRATA:
                pillars = strata_members[species][stratum]
                for scope in SCOPES:
                    if scope == "shared":
                        fg_pillars = pillars & shared_pillars
                        bg_pillars = shared_pillars
                    else:
                        fg_pillars = pillars & any_pillars
                        bg_pillars = any_pillars
                    fg_genes: set[str] = set()
                    for gi in fg_pillars:
                        for copy, genes in copy_index[species][gi].items():
                            if scope == "shared" or copy in ("A", "B"):
                                fg_genes.update(genes)
                    bg_genes_chr = set()
                    for gi in bg_pillars:
                        for genes in copy_index[species][gi].values():
                            bg_genes_chr.update(genes)
                    backgrounds["chr37_anchored"] = bg_genes_chr & genome_bg
                    fg_go = fg_genes & genome_bg
                    summary_rows.append(
                        {
                            "evidence_mode": mode,
                            "species": species,
                            "stratum": stratum,
                            "scope": scope,
                            "foreground_copies": len(fg_genes),
                            "foreground_with_GO": len(fg_go),
                            "chr37_anchored_with_GO": len(
                                backgrounds["chr37_anchored"]
                            ),
                            "genome_with_GO": len(genome_bg),
                        }
                    )
                    for bg_name, bg in backgrounds.items():
                        units = {g: go_by_species[species].get(g, set()) for g in bg}
                        index: dict[str, set] = defaultdict(set)
                        for g, terms in units.items():
                            for t in terms:
                                index[t].add(g)
                        n, N = len(fg_go), len(bg)
                        if n == 0 or N == 0:
                            continue
                        for go_id, ids in index.items():
                            k = len(fg_go & ids)
                            K = len(ids)
                            aspect, description = meta.get(go_id, ("", ""))
                            ora_rows.append(
                                {
                                    "evidence_mode": mode,
                                    "species": species,
                                    "stratum": stratum,
                                    "scope": scope,
                                    "background": bg_name,
                                    "go_id": go_id,
                                    "go_aspect": aspect,
                                    "go_description": description,
                                    "foreground_genes": n,
                                    "background_genes": N,
                                    "term_foreground": k,
                                    "term_background": K,
                                    "fold_enrichment": (k / n) / (K / N),
                                    "p_value": hypergeom_sf(k, N, K, n),
                                }
                            )
                    key = (species, stratum, scope)
                    per_species_results[key] = fg_go
                    if mode == "all":
                        for gi in sorted(fg_pillars):
                            for copy, genes in copy_index[species][gi].items():
                                for gene in sorted(genes & genome_bg):
                                    member_rows.append(
                                        {
                                            "species": species,
                                            "stratum": stratum,
                                            "scope": scope,
                                            "gtr_global_index": gi,
                                            "copy": copy,
                                            "gene_id": gene,
                                        }
                                    )

    ora = pd.DataFrame(ora_rows)
    ora["fdr"] = ora.groupby(
        ["evidence_mode", "species", "stratum", "scope", "background", "go_aspect"],
        group_keys=False,
    )["p_value"].transform(bh_group)
    ora.sort_values(
        ["evidence_mode", "stratum", "scope", "background", "species", "fdr"]
    ).to_csv(OUT / "GO_chr3_chr7_concordant_ORA.tsv", sep="\t", index=False)

    meta_rows = []
    for mode in EVIDENCE_MODES:
        for stratum in STRATA:
            for scope in SCOPES:
                for bg_name in ("genome", "chr37_anchored"):
                    sub = ora[
                        ora["evidence_mode"].eq(mode)
                        & ora["stratum"].eq(stratum)
                        & ora["scope"].eq(scope)
                        & ora["background"].eq(bg_name)
                        & ora["term_background"].ge(MIN_TERM_BG)
                    ]
                    for (go_id, aspect, desc), grp in sub.groupby(
                        ["go_id", "go_aspect", "go_description"]
                    ):
                        if grp["species"].nunique() < MIN_SPECIES:
                            continue
                        meta_rows.append(
                            {
                                "evidence_mode": mode,
                                "stratum": stratum,
                                "scope": scope,
                                "background": bg_name,
                                "go_id": go_id,
                                "go_aspect": aspect,
                                "go_description": desc,
                                "n_species_tested": grp["species"].nunique(),
                                "n_species_fdr05": int(grp["fdr"].lt(0.05).sum()),
                                "species": ",".join(sorted(grp["species"])),
                                "mean_fold": float(grp["fold_enrichment"].mean()),
                                "min_fdr": float(grp["fdr"].min()),
                                "combined_p": partial_conjunction(
                                    grp["p_value"].tolist()
                                ),
                            }
                        )
    meta_df = pd.DataFrame(meta_rows)
    if not meta_df.empty:
        meta_df["fdr"] = meta_df.groupby(
            ["evidence_mode", "stratum", "scope", "background"]
        )["combined_p"].transform(bh_group)
        meta_df.sort_values(
            ["evidence_mode", "stratum", "scope", "background", "combined_p"]
        ).to_csv(OUT / "GO_chr3_chr7_concordant_meta.tsv", sep="\t", index=False)

    pd.DataFrame(member_rows).to_csv(
        OUT / "GO_chr3_chr7_concordant_gene_members.tsv", sep="\t", index=False
    )
    pd.DataFrame(summary_rows).to_csv(
        OUT / "GO_chr3_chr7_concordant_summary.tsv", sep="\t", index=False
    )

    print("== bin classes per species (Chr03+Chr07) ==")
    print(
        trends.pivot_table(
            index="species", columns="trend_class", values="bin", aggfunc="count"
        )
        .fillna(0)
        .astype(int)
        .to_string()
    )
    print()
    print("== foreground gene copies (evidence: all) ==")
    s = pd.DataFrame(summary_rows)
    print(
        s[s["evidence_mode"].eq("all")][
            [
                "species",
                "stratum",
                "scope",
                "foreground_with_GO",
                "chr37_anchored_with_GO",
                "genome_with_GO",
            ]
        ].to_string(index=False)
    )
    print()
    if not meta_df.empty:
        sig_all = meta_df[meta_df["evidence_mode"].eq("all")]
        for (strat, scope, bg), grp in sig_all.groupby(
            ["stratum", "scope", "background"]
        ):
            sig = grp[grp["fdr"].lt(0.05)].sort_values("combined_p")
            print(
                f"== {strat} | {scope} | {bg} | terms={len(grp)} | FDR<0.05={len(sig)} "
                f"| supported by >=2 species={int(sig['n_species_fdr05'].ge(2).sum())} =="
            )
            top = sig[sig["n_species_fdr05"].ge(2)].head(12)
            if top.empty:
                top = sig.head(8)
            if top.empty:
                print("(none)")
            else:
                print(
                    top[
                        [
                            "go_description",
                            "n_species_tested",
                            "n_species_fdr05",
                            "mean_fold",
                            "combined_p",
                            "fdr",
                        ]
                    ].to_string(index=False)
                )
            print()
    else:
        print("(no term tested in >= 2 species)")


if __name__ == "__main__":
    main()
