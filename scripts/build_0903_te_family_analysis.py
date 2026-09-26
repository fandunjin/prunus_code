#!/usr/bin/env python3
"""Build the family-level TE expansion evidence package for the 0903 analysis.

The script deliberately separates:
  * unified pan-TE abundance statistics;
  * Gypsy lineage evidence from REXdb/TEsorter domains;
  * EDTA-class Mutator-like candidates;
  * orthosite, tandem-repeat and insertion-age evidence.

No subfamily names are assigned here.  Subfamily inference is prepared as a
sequence-level follow-up because the frozen pan-TE FASTA is not part of the
local working set.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PAN_ROOT = ROOT / "analysis" / "prunus_pan_te"
OUT_ROOT = ROOT / "0903" / "TE"
TABLE_ROOT = OUT_ROOT / "family_analysis"
FIG_ROOT = TABLE_ROOT / "figures"

TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
PZI_TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2"]
PSP_TARGETS = ["pspA1", "pspA2", "pspB1", "pspB2"]
PZI_A = ["pziA1", "pziA2"]
PZI_B = ["pziB1", "pziB2"]
PSP_A = ["pspA1", "pspA2"]
PSP_B = ["pspB1", "pspB2"]

SOURCES = {
    "family_summary": PAN_ROOT / "summary" / "prunus_pan_te_family_summary.tsv",
    "driver_families": PAN_ROOT / "postprocess_v2" / "te_driver_families.tsv",
    "component_decisions": PAN_ROOT / "component_decisions.tsv",
    "curated_drivers": PAN_ROOT / "followup" / "curated_summary" / "curated_te_driver_families.tsv",
    "intact_ltr": PAN_ROOT / "postprocess_v2" / "te_intact_ltr.tsv",
    "intact_age": PAN_ROOT / "followup_v2_local_flank5k" / "intact_ltr_age.tsv",
    "divergence_age": PAN_ROOT / "followup_v2_local_flank5k" / "te_divergence_age.tsv",
    "orthosite_sites": PAN_ROOT / "followup_orthosite_asm20_flank5k" / "combined" / "orthosite_site_evidence.tsv",
    "orthosite_family": PAN_ROOT / "followup_orthosite_asm20_flank5k" / "combined" / "orthosite_family_summary.tsv",
    "rexdb_domains": PAN_ROOT / "library" / "prunus_pan_te_v1.draft.fa.rexdb-plant.dom.faa",
    "sololtr_lineages": PAN_ROOT / "followup_sololtr_v2" / "sololtr_summary" / "sololtr_by_target_lineage.tsv",
}

REX_COVERAGE_MIN = 50.0
REX_PROBABILITY_MIN = 0.90
REX_EVALUE_MAX = 1e-5
RECENT_AGE_MYA = 20.0
LOW_TRF_MAX = 0.25


def ensure_dirs() -> None:
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    FIG_ROOT.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0.0)


def safe_mean(values) -> float:
    vals = [float(x) for x in values if pd.notna(x)]
    return float(np.mean(vals)) if vals else 0.0


def safe_median(values) -> float:
    vals = [float(x) for x in values if pd.notna(x)]
    return float(np.median(vals)) if vals else float("nan")


def family_base(header: str) -> str:
    match = re.search(r"(panTE_\d+(?:_(?:INT|LTR))?)", header)
    return match.group(1) if match else ""


def parse_kv_header(header: str) -> dict[str, str]:
    fields = {}
    for token in header.lstrip(">\n").split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key.strip()] = value.strip()
    return fields


def parse_rexdb(path: Path) -> pd.DataFrame:
    """Parse TEsorter/REXdb protein-domain FASTA headers."""
    records = []
    if not path.exists():
        return pd.DataFrame(columns=["pan_family", "clade", "gene", "coverage", "evalue", "probability", "header"])
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        kv = parse_kv_header(line)
        family = family_base(line)
        if not family:
            continue
        try:
            coverage = float(kv.get("coverage", "nan"))
        except ValueError:
            coverage = float("nan")
        try:
            evalue = float(kv.get("evalue", "nan"))
        except ValueError:
            evalue = float("nan")
        try:
            probability = float(kv.get("probability", "nan"))
        except ValueError:
            probability = float("nan")
        # TEsorter stores gene in the explicit key; use the name token as a fallback.
        gene = kv.get("gene", "")
        clade = kv.get("clade", "")
        records.append({
            "pan_family": family,
            "clade": clade,
            "gene": gene,
            "coverage": coverage,
            "evalue": evalue,
            "probability": probability,
            "header": line[1:],
        })
    return pd.DataFrame.from_records(records)


def summarize_rexdb(domains: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if domains.empty:
        empty = pd.DataFrame(columns=["pan_family", "rexdb_clade", "gypsy_lineage", "rexdb_domain_status"])
        return domains, empty
    domains = domains.copy()
    domains["pass"] = (
        (domains["coverage"] >= REX_COVERAGE_MIN)
        & (domains["probability"] >= REX_PROBABILITY_MIN)
        & (domains["evalue"] <= REX_EVALUE_MAX)
    )
    domains["clade"] = domains["clade"].replace("", "unclassified")
    rows = []
    for family, group in domains.groupby("pan_family", sort=True):
        passed = group[group["pass"]]
        clade_counts = Counter(passed["clade"])
        if clade_counts:
            best_clade, best_count = clade_counts.most_common(1)[0]
            tied = [k for k, v in clade_counts.items() if v == best_count]
            rex_clade = tied[0] if len(tied) == 1 else "mixed"
        else:
            best_clade, best_count, rex_clade = "", 0, ""
        genes = sorted(set(passed["gene"]))
        same_clade = passed[passed["clade"] == best_clade]
        same_genes = set(same_clade["gene"])
        strong_core = {"INT", "RT", "RH"}.issubset(same_genes)
        high_conf = (best_count >= 2 and rex_clade not in {"", "mixed"}) or strong_core
        rows.append({
            "pan_family": family,
            "rexdb_all_hit_count": len(group),
            "rexdb_pass_hit_count": len(passed),
            "rexdb_pass_clades": ",".join(sorted(clade_counts)),
            "rexdb_pass_genes": ",".join(genes),
            "rexdb_max_coverage": float(passed["coverage"].max()) if len(passed) else np.nan,
            "rexdb_min_evalue": float(passed["evalue"].min()) if len(passed) else np.nan,
            "rexdb_max_probability": float(passed["probability"].max()) if len(passed) else np.nan,
            "rexdb_clade": rex_clade,
            "rexdb_high_confidence": int(high_conf),
            "rexdb_domain_status": "passing" if len(passed) else "no_passing_hit",
        })
    return domains, pd.DataFrame(rows)


def pivot_numeric(df: pd.DataFrame, value: str) -> pd.DataFrame:
    result = df.pivot_table(index="pan_family", columns="target", values=value, aggfunc="sum", fill_value=0)
    for target in TARGET_ORDER:
        if target not in result.columns:
            result[target] = 0.0
    return result[TARGET_ORDER]


def build_master(summary: pd.DataFrame) -> pd.DataFrame:
    summary = summary.copy()
    summary["union_bp"] = num(summary["union_bp"])
    summary["annotation_count"] = num(summary["annotation_count"])
    summary["genome_fraction"] = num(summary["genome_fraction"])
    meta = summary.groupby("pan_family", as_index=True).agg(
        class_family=("class_family", lambda s: ";".join(sorted(set(s)))),
    )
    bp = pivot_numeric(summary, "union_bp")
    count = pivot_numeric(summary, "annotation_count")
    frac = pivot_numeric(summary, "genome_fraction")
    out = meta.join(bp.add_suffix("_bp")).join(count.add_suffix("_copies")).join(frac.add_suffix("_fraction"))

    def mean_cols(cols):
        return out[[f"{x}_bp" for x in cols]].mean(axis=1)

    for name, cols in {
        "pzi_a": PZI_A, "pzi_b": PZI_B, "psp_a": PSP_A, "psp_b": PSP_B,
        "pzi": PZI_TARGETS, "psp": PSP_TARGETS,
    }.items():
        out[f"{name}_mean_bp"] = mean_cols(cols)
        out[f"{name}_mean_copies"] = out[[f"{x}_copies" for x in cols]].mean(axis=1)
        out[f"{name}_mean_fraction"] = out[[f"{x}_fraction" for x in cols]].mean(axis=1)
    out["delta_a_bp"] = out["pzi_a_mean_bp"] - out["psp_a_mean_bp"]
    out["delta_b_bp"] = out["pzi_b_mean_bp"] - out["psp_b_mean_bp"]
    out["pzi_psp_excess_bp"] = out["pzi_mean_bp"] - out["psp_mean_bp"]
    out["pzi_psp_excess_fraction"] = out["pzi_mean_fraction"] - out["psp_mean_fraction"]
    out["pzi_positive_vs_matched_psp_count"] = [
        int(sum(row[f"{p}_bp"] > row[f"{s}_bp"] for p, s in zip(PZI_TARGETS, PSP_TARGETS)))
        for _, row in out.iterrows()
    ]
    pzi_values = out[[f"{x}_bp" for x in PZI_TARGETS]]
    out["pzi_haplotype_min_bp"] = pzi_values.min(axis=1)
    out["pzi_haplotype_max_bp"] = pzi_values.max(axis=1)
    out["pzi_haplotype_cv"] = pzi_values.std(axis=1, ddof=0) / pzi_values.mean(axis=1).replace(0, np.nan)
    out["pzi_haplotype_cv"] = out["pzi_haplotype_cv"].fillna(0)
    out["pzi_all_haplotypes_present"] = (pzi_values > 0).all(axis=1).astype(int)
    out["pzi_min_max_ratio"] = out["pzi_haplotype_min_bp"] / out["pzi_haplotype_max_bp"].replace(0, np.nan)
    out["pzi_min_max_ratio"] = out["pzi_min_max_ratio"].fillna(0)
    out["pzi_consistency"] = np.select(
        [
            (out["pzi_positive_vs_matched_psp_count"] == 4) & (out["pzi_min_max_ratio"] >= 0.5),
            (out["pzi_all_haplotypes_present"] == 1) & (out["pzi_min_max_ratio"] >= 0.25),
        ],
        ["consistent_expansion", "partial_consistency"],
        default="inconsistent_or_absent",
    )
    out = out.reset_index()
    return out


def aggregate_optional_evidence(master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Join curated, LTR age, orthosite and manual TRF evidence."""
    out = master.copy()

    decisions = read_tsv(SOURCES["component_decisions"], required=False)
    curated = read_tsv(SOURCES["curated_drivers"], required=False)
    if not decisions.empty:
        keep = [c for c in ["pan_family", "positive_rank", "class_family", "positive_excess_bp", "pzi_mean_tandem_overlap_fraction", "automated_status", "manual_decision", "curated_family", "manual_evidence", "notes"] if c in decisions]
        decisions = decisions[keep].copy()
        decisions = decisions.rename(columns={"class_family": "audited_class_family", "manual_decision": "audit_manual_decision", "curated_family": "audit_curated_family", "manual_evidence": "audit_manual_evidence"})
        out = out.merge(decisions, on="pan_family", how="left")
    else:
        out["pzi_mean_tandem_overlap_fraction"] = np.nan

    if not curated.empty:
        curated_keep = [c for c in ["pan_family", "manual_decision", "manual_evidence", "audit_notes", "curated_rank", "pzi_intact_ltr_count", "psp_intact_ltr_count", "pzi_intact_ltr_consistent", "psp_intact_ltr_consistent"] if c in curated]
        curated = curated[curated_keep].copy().rename(columns={"manual_decision": "curated_manual_decision", "manual_evidence": "curated_manual_evidence"})
        out = out.merge(curated, on="pan_family", how="left")
    for col in ["pzi_intact_ltr_count", "psp_intact_ltr_count", "pzi_intact_ltr_consistent", "psp_intact_ltr_consistent"]:
        if col not in out:
            out[col] = 0
        out[col] = num(out[col])
    if "pzi_mean_tandem_overlap_fraction" not in out:
        out["pzi_mean_tandem_overlap_fraction"] = np.nan
    out["pzi_mean_tandem_overlap_fraction"] = pd.to_numeric(out["pzi_mean_tandem_overlap_fraction"], errors="coerce")
    out["manual_decision_final"] = out.get("curated_manual_decision", pd.Series("", index=out.index)).replace("", np.nan).fillna(out.get("audit_manual_decision", pd.Series("", index=out.index)))
    out["manual_evidence_final"] = out.get("curated_manual_evidence", pd.Series("", index=out.index)).replace("", np.nan).fillna(out.get("audit_manual_evidence", pd.Series("", index=out.index)))

    intact = read_tsv(SOURCES["intact_age"], required=False)
    if intact.empty:
        intact = read_tsv(SOURCES["intact_ltr"], required=False)
    if not intact.empty and "pan_family" in intact.columns:
        intact["estimated_age_mya"] = pd.to_numeric(intact.get("estimated_age_mya", np.nan), errors="coerce")
        age_rows = []
        for family, g in intact.groupby("pan_family"):
            pzi = g[g.get("group", "") == "PZI"] if "group" in g else g[g["target"].isin(PZI_TARGETS)]
            psp = g[g.get("group", "") == "PSP"] if "group" in g else g[g["target"].isin(PSP_TARGETS)]
            ages = pzi["estimated_age_mya"].dropna()
            age_rows.append({
                "pan_family": family,
                "pzi_intact_age_record_count": len(pzi),
                "psp_intact_age_record_count": len(psp),
                "pzi_intact_age_median_mya": float(ages.median()) if len(ages) else np.nan,
                "pzi_intact_recent_age_fraction": float((ages <= RECENT_AGE_MYA).mean()) if len(ages) else np.nan,
            })
        out = out.merge(pd.DataFrame(age_rows), on="pan_family", how="left")
    else:
        for col in ["pzi_intact_age_record_count", "psp_intact_age_record_count", "pzi_intact_age_median_mya", "pzi_intact_recent_age_fraction"]:
            out[col] = np.nan

    div = read_tsv(SOURCES["divergence_age"], required=False)
    if not div.empty and "pan_family" in div.columns:
        div["estimated_age_mya"] = pd.to_numeric(div.get("estimated_age_mya", np.nan), errors="coerce")
        div["exclusive_bp"] = num(div.get("exclusive_bp", pd.Series(dtype=float)))
        div_rows = []
        for family, g in div.groupby("pan_family"):
            pzi = g[g.get("group", "") == "PZI"] if "group" in g else g
            ages = pzi["estimated_age_mya"].dropna()
            weights = pzi.loc[ages.index, "exclusive_bp"] if len(ages) else pd.Series(dtype=float)
            weighted = np.average(ages, weights=np.maximum(weights, 1)) if len(ages) else np.nan
            div_rows.append({
                "pan_family": family,
                "pzi_divergence_age_median_mya": float(ages.median()) if len(ages) else np.nan,
                "pzi_divergence_age_weighted_mya": float(weighted) if len(ages) else np.nan,
                "pzi_divergence_recent_fraction": float((ages <= RECENT_AGE_MYA).mean()) if len(ages) else np.nan,
            })
        out = out.merge(pd.DataFrame(div_rows), on="pan_family", how="left")
    else:
        for col in ["pzi_divergence_age_median_mya", "pzi_divergence_age_weighted_mya", "pzi_divergence_recent_fraction"]:
            out[col] = np.nan

    ortho = read_tsv(SOURCES["orthosite_sites"], required=False)
    if not ortho.empty and "pan_family" in ortho.columns:
        rows = []
        for family, g in ortho.groupby("pan_family"):
            ec = g.get("evidence_class", pd.Series(dtype=str)).astype(str)
            rows.append({
                "pan_family": family,
                "orthosite_site_count": len(g),
                "orthosite_resolved_count": int((ec != "unresolved").sum()),
                "orthosite_pzi_specific_count": int(ec.str.contains("pzi_specific", na=False).sum()),
                "orthosite_shared_count": int(ec.str.contains("shared_pzi_psp", na=False).sum()),
                "orthosite_status": "supported" if (ec != "unresolved").any() else "unresolved_only",
            })
        out = out.merge(pd.DataFrame(rows), on="pan_family", how="left")
    else:
        for col in ["orthosite_site_count", "orthosite_resolved_count", "orthosite_pzi_specific_count", "orthosite_shared_count"]:
            out[col] = 0
        out["orthosite_status"] = "not_available"
    out[["orthosite_site_count", "orthosite_resolved_count", "orthosite_pzi_specific_count", "orthosite_shared_count"]] = out[["orthosite_site_count", "orthosite_resolved_count", "orthosite_pzi_specific_count", "orthosite_shared_count"]].fillna(0)
    out["orthosite_status"] = out["orthosite_status"].fillna("not_available")

    out["pzi_age_median_mya"] = out["pzi_intact_age_median_mya"].fillna(out["pzi_divergence_age_median_mya"])
    out["pzi_recent_age_fraction"] = out["pzi_intact_recent_age_fraction"].fillna(out["pzi_divergence_recent_fraction"])
    out["low_trf_overlap"] = np.where(
        out["pzi_mean_tandem_overlap_fraction"].notna(),
        out["pzi_mean_tandem_overlap_fraction"] <= LOW_TRF_MAX,
        np.nan,
    )
    out["intact_or_orthosite_support"] = (
        (out["pzi_intact_ltr_count"] > 0) | (out["orthosite_resolved_count"] > 0)
    ).astype(int)
    out["recent_insertion_support"] = (out["pzi_recent_age_fraction"] >= 0.5).astype(int)
    out["evidence_support_score"] = (
        (out["pzi_consistency"] == "consistent_expansion").astype(int)
        + out["low_trf_overlap"].fillna(False).astype(int)
        + out["intact_or_orthosite_support"]
        + out["recent_insertion_support"]
    )
    out["evidence_tier"] = np.select(
        [
            (out["evidence_support_score"] >= 4) & (out["manual_decision_final"].fillna("") != "exclude"),
            (out["evidence_support_score"] >= 2) & (out["manual_decision_final"].fillna("") != "exclude"),
        ],
        ["high_support", "moderate_support"],
        default="exploratory_or_unresolved",
    )
    return out, decisions, curated


def add_family_classifications(master: pd.DataFrame, rex_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = master.merge(rex_summary, on="pan_family", how="left")
    for col in ["rexdb_clade", "rexdb_pass_clades", "rexdb_pass_genes", "rexdb_domain_status"]:
        if col not in out:
            out[col] = ""
        out[col] = out[col].fillna("")
    out["rexdb_high_confidence"] = num(out.get("rexdb_high_confidence", pd.Series(0, index=out.index)))
    # A strong REXdb clade can refine an EDTA LTR/unknown family.  An EDTA
    # Gypsy family without passing domains remains Gypsy/unclassified.
    is_ltr = out["class_family"].str.startswith("LTR/")
    has_strong_gypsy = (out["rexdb_high_confidence"] == 1) & out["rexdb_clade"].isin(
        ["Tekay", "Athila", "CRM", "Retand", "Ogre", "Galadriel", "Reina", "GMR", "OTA"]
    )
    out["gypsy_lineage"] = np.where(
        is_ltr & has_strong_gypsy,
        "Gypsy/" + out["rexdb_clade"],
        np.where(is_ltr & out["class_family"].str.startswith("LTR/Gypsy"), "Gypsy/unclassified", "not_applicable"),
    )
    out["mutator_class"] = np.select(
        [
            out["class_family"].str.startswith(("DNA/DTM", "MITE/DTM")),
            out["class_family"].str.startswith("DNA/MULE-MuDR"),
        ],
        ["Mutator-like/DTM", "MULE-like/MuDR"],
        default="not_applicable",
    )
    out["mutator_domain_status"] = np.where(
        out["mutator_class"] == "not_applicable", "not_applicable", "not_available_in_current_REXdb_table"
    )
    out["mutator_tir_status"] = np.where(
        out["mutator_class"] == "not_applicable", "not_applicable", "not_available_in_current_tables"
    )
    out["mutator_intact_element_status"] = np.where(
        out["mutator_class"] == "not_applicable", "not_applicable", "not_available_in_current_tables"
    )
    out["mutator_claim_status"] = np.select(
        [out["mutator_class"] == "Mutator-like/DTM", out["mutator_class"] == "MULE-like/MuDR"],
        ["exploratory_DTM_only", "exploratory_MuDR_class_only"],
        default="not_applicable",
    )
    out["mutator_evidence_note"] = np.where(
        out["mutator_class"] == "Mutator-like/DTM",
        "EDTA DTM class; transposase/TIR/intact evidence must be confirmed separately",
        np.where(out["mutator_class"] == "MULE-like/MuDR", "EDTA MuDR class; domain/TIR confirmation required", ""),
    )
    mutator = out[out["mutator_class"] != "not_applicable"].copy()
    return out, mutator


def prepare_candidate_table(master: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "pan_family", "class_family", "pzi_psp_excess_bp", "pzi_psp_excess_fraction",
        "pzi_consistency", "pzi_positive_vs_matched_psp_count", "pzi_haplotype_cv",
        "pzi_mean_tandem_overlap_fraction", "low_trf_overlap", "pzi_intact_ltr_count",
        "orthosite_resolved_count", "pzi_age_median_mya", "pzi_recent_age_fraction",
        "rexdb_clade", "gypsy_lineage", "mutator_class", "evidence_support_score", "evidence_tier",
        "manual_decision_final", "manual_evidence_final", "audit_notes",
    ]
    cols = [c for c in cols if c in master.columns]
    cand = master[cols].copy()
    cand = cand.sort_values(["pzi_psp_excess_bp", "evidence_support_score"], ascending=[False, False])
    cand["candidate_rank"] = np.arange(1, len(cand) + 1)
    return cand


def save_pub(fig: plt.Figure, base: Path) -> None:
    fig.savefig(str(base.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(base.with_suffix(".pdf")), bbox_inches="tight")
    fig.savefig(str(base.with_suffix(".png")), dpi=600, bbox_inches="tight")
    fig.savefig(str(base.with_suffix(".tiff")), dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(master: pd.DataFrame) -> None:
    # The primary heatmap is an evidence plot: manually excluded tandem-array
    # artifacts stay in the master table but are not allowed to dominate it.
    d = master[master["manual_decision_final"].fillna("") != "exclude"].sort_values("pzi_psp_excess_bp", ascending=False).head(30).copy()
    d.to_csv(TABLE_ROOT / "figure_heatmap_data.tsv", sep="\t", index=False)
    mat = np.log10(d[[f"{x}_bp" for x in TARGET_ORDER]].to_numpy(dtype=float) + 1.0)
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    im = ax.imshow(mat, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_xticks(range(len(TARGET_ORDER)), [x.replace("pzi", "PZI-").replace("psp", "PSP-").replace("ppa", "PPA-").upper() for x in TARGET_ORDER], rotation=50, ha="right")
    ax.set_yticks(range(len(d)), d["pan_family"].tolist(), fontsize=5.8)
    ax.set_xlabel("Target genome")
    ax.set_ylabel("Pan-TE family (ranked by PZI–PSP excess)")
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
    cbar.set_label("log10(union bp + 1)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_pub(fig, FIG_ROOT / "te_family_expansion_heatmap")


def plot_driver_scatter(master: pd.DataFrame) -> None:
    d = master[master["pzi_mean_tandem_overlap_fraction"].notna()].copy()
    if d.empty:
        return
    palette = {"LTR/Gypsy": "#8C2D04", "LTR/unknown": "#4A1486", "DNA/DTM": "#006D77", "MITE/DTM": "#2C7FB8", "DNA/CMC-EnSpm": "#6A51A3", "DNA/Helitron": "#3F007D"}
    d.to_csv(TABLE_ROOT / "figure_driver_scatter_data.tsv", sep="\t", index=False)
    excluded = d["manual_decision_final"].fillna("") == "exclude"
    colors = ["#BDBDBD" if x else next((v for k, v in palette.items() if str(c).startswith(k)), "#444444") for x, c in zip(excluded, d["class_family"])]
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.scatter(d.loc[~excluded, "pzi_mean_tandem_overlap_fraction"], d.loc[~excluded, "pzi_psp_excess_bp"] / 1e6, c=np.asarray(colors)[~excluded.to_numpy()], s=22, alpha=0.82, edgecolors="white", linewidths=0.25, label="not excluded")
    if excluded.any():
        ax.scatter(d.loc[excluded, "pzi_mean_tandem_overlap_fraction"], d.loc[excluded, "pzi_psp_excess_bp"] / 1e6, c="#BDBDBD", s=28, alpha=0.8, marker="x", linewidths=0.7, label="manual tandem-repeat exclusion")
    selected = d[(~excluded) & (d["pzi_mean_tandem_overlap_fraction"] <= LOW_TRF_MAX)].sort_values("pzi_psp_excess_bp", ascending=False).head(8)
    for _, row in selected.iterrows():
        ax.annotate(row["pan_family"], (row["pzi_mean_tandem_overlap_fraction"], row["pzi_psp_excess_bp"] / 1e6), xytext=(3, 3), textcoords="offset points", fontsize=5.7)
    ax.axvline(LOW_TRF_MAX, color="#222222", lw=0.8, ls="--")
    ax.axhline(0, color="#777777", lw=0.6)
    ax.set_xlabel("PZI mean TRF overlap fraction")
    ax.set_ylabel("PZI–PSP excess (Mb)")
    ax.set_xlim(left=-0.02)
    ax.legend(loc="upper left", fontsize=5.5, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_pub(fig, FIG_ROOT / "te_driver_evidence_scatter")


def plot_lineage_mutator(master: pd.DataFrame) -> None:
    valid = master[master["manual_decision_final"].fillna("") != "exclude"]
    gy = valid[valid["gypsy_lineage"].str.startswith("Gypsy/")].copy()
    gy = gy[gy["pzi_psp_excess_bp"] > 0]
    gy = gy.groupby("gypsy_lineage", as_index=False)["pzi_psp_excess_bp"].sum().sort_values("pzi_psp_excess_bp", ascending=True)
    mut = valid[valid["mutator_class"] != "not_applicable"].copy().sort_values("pzi_psp_excess_bp", ascending=True).tail(12)
    pd.concat([
        gy.assign(plot_group="Gypsy lineage"),
        mut[["pan_family", "mutator_class", "pzi_psp_excess_bp"]].assign(plot_group="DTM/MuDR candidate"),
    ], ignore_index=True, sort=False).to_csv(TABLE_ROOT / "figure_lineage_mutator_data.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), gridspec_kw={"width_ratios": [1.0, 1.45]})
    if not gy.empty:
        axes[0].barh(gy["gypsy_lineage"], gy["pzi_psp_excess_bp"] / 1e6, color="#7F0000")
    axes[0].set_xlabel("PZI–PSP excess (Mb)")
    axes[0].set_ylabel("Gypsy lineage")
    axes[0].spines[["top", "right"]].set_visible(False)
    if not mut.empty:
        colors = ["#006D77" if x == "Mutator-like/DTM" else "#1D3557" for x in mut["mutator_class"]]
        axes[1].barh(mut["pan_family"], mut["pzi_psp_excess_bp"] / 1e6, color=colors)
    axes[1].set_xlabel("PZI–PSP excess (Mb)")
    axes[1].set_ylabel("DTM / MuDR candidate family")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout(w_pad=2.2)
    save_pub(fig, FIG_ROOT / "te_gypsy_mutator_summary")


def write_lineage_summaries(master: pd.DataFrame) -> None:
    valid = master[master["manual_decision_final"].fillna("") != "exclude"].copy()
    gy = valid[valid["gypsy_lineage"].str.startswith("Gypsy/")]
    if gy.empty:
        gy_summary = pd.DataFrame(columns=["gypsy_lineage", "family_count", "pzi_psp_excess_bp", "consistent_family_count", "high_support_family_count"])
    else:
        gy_summary = gy.groupby("gypsy_lineage", as_index=False).agg(
            family_count=("pan_family", "nunique"),
            pzi_psp_excess_bp=("pzi_psp_excess_bp", "sum"),
            consistent_family_count=("pzi_consistency", lambda x: int((x == "consistent_expansion").sum())),
            high_support_family_count=("evidence_tier", lambda x: int((x == "high_support").sum())),
        ).sort_values("pzi_psp_excess_bp", ascending=False)
    gy_summary.to_csv(TABLE_ROOT / "gypsy_lineage_summary.tsv", sep="\t", index=False)
    mut = valid[valid["mutator_class"] != "not_applicable"]
    if mut.empty:
        mut_summary = pd.DataFrame(columns=["mutator_class", "family_count", "pzi_psp_excess_bp", "consistent_family_count", "low_trf_audited_count"])
    else:
        mut_summary = mut.groupby("mutator_class", as_index=False).agg(
            family_count=("pan_family", "nunique"),
            pzi_psp_excess_bp=("pzi_psp_excess_bp", "sum"),
            consistent_family_count=("pzi_consistency", lambda x: int((x == "consistent_expansion").sum())),
            low_trf_audited_count=("low_trf_overlap", lambda x: int((x == True).sum())),
        ).sort_values("pzi_psp_excess_bp", ascending=False)
    mut_summary.to_csv(TABLE_ROOT / "mutator_class_summary.tsv", sep="\t", index=False)


def write_subfamily_readme(master: pd.DataFrame) -> None:
    eligible = master[
        ((master["gypsy_lineage"].str.startswith("Gypsy/") & (master["pzi_psp_excess_bp"] > 0))
        | (master["mutator_class"] != "not_applicable"))
    ].copy()
    eligible = eligible[["pan_family", "class_family", "gypsy_lineage", "mutator_class", "pzi_psp_excess_bp", "evidence_tier"]].sort_values("pzi_psp_excess_bp", ascending=False)
    eligible.to_csv(TABLE_ROOT / "subfamily_analysis_manifest.tsv", sep="\t", index=False)
    valid = master[master["manual_decision_final"].fillna("") != "exclude"].copy()
    gypsy_tree = valid[
        valid["gypsy_lineage"].str.startswith("Gypsy/")
        & (valid["rexdb_high_confidence"] == 1)
        & (valid["pzi_psp_excess_bp"] > 0)
        & valid["evidence_tier"].isin(["high_support", "moderate_support"])
    ].copy()
    gypsy_tree["tree_mode"] = "Gypsy_RT_RH_INT"
    gypsy_tree["selection_reason"] = "named REXdb lineage; positive PZI-PSP excess; non-excluded moderate/high evidence"
    mutator_tree = valid[
        (valid["mutator_class"] != "not_applicable")
        & (valid["pzi_psp_excess_bp"] > 0)
        & (valid["manual_decision_final"].fillna("") == "include")
    ].copy()
    if len(mutator_tree) < 2:
        fallback = valid[
            (valid["mutator_class"] != "not_applicable")
            & (valid["pzi_psp_excess_bp"] > 0)
            & valid["evidence_tier"].isin(["high_support", "moderate_support"])
        ].copy().sort_values("pzi_psp_excess_bp", ascending=False)
        mutator_tree = pd.concat([mutator_tree, fallback], ignore_index=True).drop_duplicates("pan_family")
    mutator_tree = mutator_tree.sort_values("pzi_psp_excess_bp", ascending=False).head(10)
    mutator_tree["tree_mode"] = "Mutator_transposase"
    mutator_tree["selection_reason"] = "DTM/MuDR class; positive PZI-PSP excess; transposase/TIR validation required"
    tree_selection = pd.concat([gypsy_tree, mutator_tree], ignore_index=True, sort=False)
    tree_selection = tree_selection[[
        "pan_family", "class_family", "gypsy_lineage", "mutator_class", "tree_mode",
        "pzi_psp_excess_bp", "pzi_consistency", "pzi_mean_tandem_overlap_fraction",
        "evidence_tier", "manual_decision_final", "selection_reason",
    ]].sort_values(["tree_mode", "pzi_psp_excess_bp"], ascending=[True, False])
    tree_selection.to_csv(TABLE_ROOT / "subfamily_tree_selection.tsv", sep="\t", index=False)
    text = f"""# Subfamily analysis preparation

This manifest contains {len(eligible)} candidate pan-families. `subfamily_tree_selection.tsv` is the executable first-pass set: {len(gypsy_tree)} named Gypsy families and {len(mutator_tree)} DTM/MuDR candidates. It is intentionally a preparation package, not a subfamily assignment. The local working directory does not contain the frozen pan-TE FASTA or the 12 target mask genomes required to extract family-level copies.

## HPC execution

Copy `subfamily_tree_selection.tsv` to the HPC project, then run:

```bash
THREADS=24 jsub < jobs/haploid_30_prunus_te_subfamily_trees.sh
```

The job verifies the SHA256 recorded in `frozen_library_manifest.tsv`, extracts copy-level intervals from all 12 standardized pan-library RepeatMasker `.out` files, uses `bedtools getfasta -s`, performs MMseqs2 or CD-HIT clustering, runs TEsorter on Gypsy copy FASTAs, and sends validated RT/RH/INT domain FASTAs to MAFFT and IQ-TREE. It writes `subfamily_tree_status.tsv`, `tree_supported_clades.tsv` and `subfamily_calls.tsv`; a tree is marked complete only if the required domain sequence set and both IQ-TREE support metrics are present.

## Required input

Use the frozen `panTE` FASTA used for the final unified mask, or export one sequence FASTA per `pan_family`. Keep the original family ID in every FASTA header.

## Gypsy candidates

1. Filter families with `gypsy_lineage` support in `gypsy_rexdb_evidence.tsv` and retain only passing REXdb hits (coverage >= {REX_COVERAGE_MIN:g}, probability >= {REX_PROBABILITY_MIN:g}, E-value <= {REX_EVALUE_MAX:g}).
2. For each family, cluster full-length copies with MMseqs2 or CD-HIT.
3. Extract RT/RH (and INT where available) from the TEsorter/REXdb protein-domain output. Do not use the family consensus alone as a phylogenetic sequence.
4. Align each domain with MAFFT and infer a tree with IQ-TREE.

Example commands after the frozen FASTA and domain table are available:

```bash
mmseqs easy-cluster FAMILY.fa FAMILY_mmseqs FAMILY_tmp --min-seq-id 0.80 -c 0.80 --cov-mode 1
mafft --auto FAMILY_RT_RH.faa > FAMILY_RT_RH.aln.faa
iqtree2 -s FAMILY_RT_RH.aln.faa -m MFP -B 1000 --alrt 1000 -T AUTO --prefix FAMILY_RT_RH
```

Define a subfamily only when a clade has UFboot >= 95 and SH-aLRT >= 80, contains at least three independent copies, and is supported by the same lineage/domain assignment. Report the number of copies, alignment length, and both support values.

## Mutator-like candidates

Use EDTA `DNA/DTM` or `MITE/DTM` as the initial class only. The job stops a Mutator tree at `missing_mutator_domain` unless an EDTA/HMMER transposase FASTA is supplied with `MUTATOR_DOMAIN_ROOT`. Confirm a Mutator/MULE transposase domain, terminal inverted repeats and, where possible, an intact element before using a subfamily name. Do not apply Gypsy REXdb logic to DTM candidates. Extract the transposase region, align with MAFFT and infer an IQ-TREE tree using the same support thresholds.

## Interpretation rule

Families without domain or terminal-structure support remain `Gypsy/unclassified` or `Mutator-like/DTM exploratory`; they must not be given a named subfamily.
"""
    (TABLE_ROOT / "subfamily_analysis_README.md").write_text(text, encoding="utf-8")


def write_results_summary(master: pd.DataFrame) -> None:
    total_path = PAN_ROOT / "summary" / "prunus_pan_te_total_summary.tsv"
    total = read_tsv(total_path, required=False)
    lines = [
        "# Unified pan-TE family analysis",
        "",
        "This package separates abundance from biological validation. All family-level abundance values come from the frozen unified pan-TE family summary. Manually excluded high-TRF families remain in the master table but are excluded from the primary evidence figures.",
        "",
        "## Species-level context",
        "",
    ]
    if not total.empty:
        total["genome_bp"] = num(total["genome_bp"])
        total["repeat_union_bp"] = num(total["repeat_union_bp"])
        total["repeat_union_fraction"] = num(total["repeat_union_fraction"])
        if "repeatmasker_out" in total:
            total["repeatmasker_out"] = total["repeatmasker_out"].astype(str).str.replace("\\\\", "/").str.split("/").str[-1]
        total.to_csv(TABLE_ROOT / "species_TE_context.tsv", sep="\t", index=False)
        for group, g in total.groupby("group"):
            lines.append(
                f"- {group}: mean genome {g['genome_bp'].mean()/1e6:.1f} Mb; mean unified TE union {g['repeat_union_bp'].mean()/1e6:.1f} Mb; mean TE fraction {g['repeat_union_fraction'].mean():.3f}."
            )
    lines += [
        "",
        "## Family evidence rules",
        "",
        f"- PZI consistency: all four PZI haplotypes exceed their matched PSP haplotypes and the PZI family abundance min/max ratio is at least 0.5.",
        f"- Low tandem-repeat overlap: audited PZI mean TRF overlap <= {LOW_TRF_MAX:.2f}; unaudited families are not treated as low-overlap.",
        f"- Recent insertion support: at least half of the available PZI age records are <= {RECENT_AGE_MYA:.0f} Ma.",
        f"- Gypsy lineage support: REXdb coverage >= {REX_COVERAGE_MIN:.0f}, probability >= {REX_PROBABILITY_MIN:.2f}, E-value <= {REX_EVALUE_MAX:g}; at least two concordant passing domains or an INT+RT+RH combination. Otherwise the family is Gypsy/unclassified.",
        "- Mutator-like families use EDTA DTM or MuDR classes only. REXdb Gypsy logic is not applied to DTM families; transposase, TIR and intact-element confirmation remains a separate requirement.",
        "",
        "## Highest-priority retained families",
        "",
    ]
    valid = master[master["manual_decision_final"].fillna("") != "exclude"].sort_values("pzi_psp_excess_bp", ascending=False).head(12)
    for _, row in valid.iterrows():
        lineage = row.get("gypsy_lineage", "not_applicable")
        mclass = row.get("mutator_class", "not_applicable")
        lines.append(
            f"- `{row['pan_family']}` ({row['class_family']}): PZI-PSP excess {row['pzi_psp_excess_bp']/1e6:.2f} Mb; consistency `{row['pzi_consistency']}`; evidence tier `{row['evidence_tier']}`; Gypsy `{lineage}`; Mutator `{mclass}`."
        )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "The current package supports family-level expansion and evidence prioritization. It does not assign Gypsy or Mutator subfamily names. Subfamily claims require copy-level clustering, domain extraction, alignment and supported IQ-TREE clades using the frozen pan-TE FASTA; see `subfamily_analysis_README.md`.",
    ]
    (TABLE_ROOT / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    summary = read_tsv(SOURCES["family_summary"])
    master = build_master(summary)
    domains = parse_rexdb(SOURCES["rexdb_domains"])
    domains, rex_summary = summarize_rexdb(domains)
    master, _, _ = aggregate_optional_evidence(master)
    master, mutator = add_family_classifications(master, rex_summary)
    master = master.sort_values("pzi_psp_excess_bp", ascending=False)

    # Keep a compact but complete family-level table plus the raw target columns.
    master.to_csv(TABLE_ROOT / "pan_family_master.tsv", sep="\t", index=False, na_rep="NA")
    domains.to_csv(TABLE_ROOT / "gypsy_rexdb_domain_hits.tsv", sep="\t", index=False, na_rep="NA")
    gypsy = master[master["gypsy_lineage"].str.startswith("Gypsy/")].copy()
    gypsy.to_csv(TABLE_ROOT / "gypsy_rexdb_evidence.tsv", sep="\t", index=False, na_rep="NA")
    mutator.to_csv(TABLE_ROOT / "mutator_candidate_evidence.tsv", sep="\t", index=False, na_rep="NA")
    prepare_candidate_table(master).to_csv(TABLE_ROOT / "candidate_family_evidence.tsv", sep="\t", index=False, na_rep="NA")

    plot_heatmap(master)
    plot_driver_scatter(master)
    plot_lineage_mutator(master)
    write_lineage_summaries(master)
    write_subfamily_readme(master)
    write_results_summary(master)

    manifest = {
        "output_root": "0903/TE",
        "target_order": TARGET_ORDER,
        "source_files": {k: str(v.relative_to(ROOT)) for k, v in SOURCES.items()},
        "rexdb_thresholds": {
            "coverage_min": REX_COVERAGE_MIN,
            "probability_min": REX_PROBABILITY_MIN,
            "evalue_max": REX_EVALUE_MAX,
        },
        "evidence_rules": {
            "low_trf_overlap_max": LOW_TRF_MAX,
            "recent_age_mya": RECENT_AGE_MYA,
            "pzi_consistent_expansion": "all four PZI haplotypes exceed matched PSP haplotypes and min/max PZI family abundance >= 0.5",
            "gypsy_without_domain": "Gypsy/unclassified",
            "mutator": "EDTA DTM/MuDR class; no Gypsy REXdb logic applied",
        },
        "outputs": [
            "family_analysis/pan_family_master.tsv",
            "family_analysis/gypsy_rexdb_domain_hits.tsv",
            "family_analysis/gypsy_rexdb_evidence.tsv",
            "family_analysis/mutator_candidate_evidence.tsv",
            "family_analysis/gypsy_lineage_summary.tsv",
            "family_analysis/mutator_class_summary.tsv",
            "family_analysis/candidate_family_evidence.tsv",
            "family_analysis/subfamily_analysis_manifest.tsv",
            "family_analysis/subfamily_tree_selection.tsv",
            "family_analysis/subfamily_analysis_README.md",
            "family_analysis/analysis_summary.md",
            "family_analysis/species_TE_context.tsv",
            "family_analysis/figures/te_family_expansion_heatmap.{svg,pdf,png,tiff}",
            "family_analysis/figures/te_driver_evidence_scatter.{svg,pdf,png,tiff}",
            "family_analysis/figures/te_gypsy_mutator_summary.{svg,pdf,png,tiff}",
        ],
        "note": "Subfamily names are not assigned until sequence-level clustering, domain extraction, alignment and phylogeny are completed with the frozen pan-TE FASTA.",
    }
    (TABLE_ROOT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote family analysis to {TABLE_ROOT}")
    print(f"Families: {len(master)}; Gypsy families: {len(gypsy)}; DTM/MuDR candidates: {len(mutator)}")


if __name__ == "__main__":
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    })
    main()
