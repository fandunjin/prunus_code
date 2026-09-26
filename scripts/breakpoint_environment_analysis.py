#!/usr/bin/env python3
"""Breakpoint gene-environment analysis for the racemose A/B rearrangement story.

For every WGDI junction candidate in PZI/PSP/PPA, classify the window by
subgenome track (A/B) and by whether the same source-unit junction is also
detected in PCAM and/or PPE (conserved vs lineage-restricted), then test:

  * gene content and inflorescence/meristem candidate genes vs
    chromosome-matched random windows (empirical permutation);
  * TE content (te, gypsy, copia, helitron, tir coverage) vs the same null;
  * family-level breakdown for the candidate gene families.

Outputs land in analysis/breakpoint_environment_20260912/.
"""
from __future__ import annotations

import csv
import importlib.util
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "0903" / "karyotype" / "fig4_breakpoint_candidates"
OUT = ROOT / "analysis" / "breakpoint_environment_20260912"
RNG = np.random.default_rng(20260912)
N_RANDOM = 1000
WINDOWS = {"main": 200_000, "narrow": 50_000, "wide": 500_000}
TE_TRACKS = ["te_coverage", "gypsy_coverage", "copia_coverage", "helitron_coverage", "tir_coverage", "gene_density"]

SPECIES_LABELS = {
    "PZI": ["pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2"],
    "PSP": ["psp_hapA1", "psp_hapA2", "psp_hapB1", "psp_hapB2"],
    "PPA": ["ppa_hapA", "ppa_hapB"],
}
# Strict Swiss-Prot description patterns (avoid the AGL*/RAP2/ERF false
# positives that keyword screens are prone to; see the Fig. 5 v1 lesson).
CANDIDATES = [
    ("LFY", r"Floricaula|LEAFY protein|\bLEAFY\b"),
    ("AP1", r"Floral homeotic protein APETALA 1|APETALA1\b"),
    ("FUL", r"FRUITFULL|\bFUL\b"),
    ("TFL1/CEN", r"Terminal flower 1\b|TERMINAL FLOWER 1\b"),
    ("FT", r"Flowering locus T|\bFT protein\b"),
    ("SOC1", r"SUPPRESSOR OF OVEREXPRESSION|MADS-box protein SOC1"),
    ("SVP", r"SHORT VEGETATIVE PHASE"),
    ("AP2", r"Floral homeotic protein APETALA 2"),
    ("SEP", r"SEPALLATA"),
    ("AGL6", r"MADS-box protein AGL6\b"),
    ("AG", r"Floral homeotic protein AGAMOUS(?!-like)"),
    ("WUS", r"WUSCHEL"),
    ("CLV", r"CLAVATA"),
    ("BRC1", r"BRANCHED 1\b|BRC1\b"),
    ("LAX", r"LAX PANICLE|\bLAX\b"),
    ("KNOX/BELL", r"knotted-1-like|Homeobox protein KNAT|KNOTTED|BELL homeodomain"),
    ("PIN", r"PIN-FORMED|PINOID"),
    ("TCP", r"TEOSINTE BRANCHED|TCP\d+"),
    ("SPL", r"SQUAMOSA PROMOTER|SPL\d+"),
]
FAMILY_RES = [(name, re.compile(pattern, re.I)) for name, pattern in CANDIDATES]
EXPANSION_FAMILIES = ["SPL9/15", "SVP", "TCP", "BRC1", "AP1"]
# SPL9/15 etc. are subfamily-level; for the keyword screen we tag the parent
# family and additionally record the SPL9/15 / TCP2-4 / DAM-like subfamily via
# description text where available.


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def load_track_module():
    spec = importlib.util.spec_from_file_location(
        "bpw", ROOT / "scripts" / "analyze_fig4_breakpoint_windows.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_index(label: str, gff: Path) -> dict[str, list[tuple[float, int, int]]]:
    """chromosome -> list of (midpoint, start, end) at gene level."""
    genes: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    with gff.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            gene = re.sub(r"-mRNA-\d+$", "", f[1])
            try:
                start, end = int(f[2]), int(f[3])
            except ValueError:
                continue
            old = genes[f[0]].get(gene)
            if old is None or (end - start) > (old[1] - old[0]):
                genes[f[0]][gene] = (start, end)
    out: dict[str, list[tuple[float, int, int]]] = {}
    for chrom, gmap in genes.items():
        rows = []
        for _gene, (start, end) in gmap.items():
            rows.append(((start + end) / 2.0, start, end))
        rows.sort()
        out[chrom] = rows
    return out


def functional_families(root: Path, species: str):
    path = root / "0903" / "annotation" / "functional" / species / "swissprot_functional_annotation.tsv"
    families: dict[str, set[str]] = {}
    descriptions: dict[str, str] = {}
    if not path.is_file():
        return families, descriptions
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ident = row.get("protein_id", "")
            if not ident:
                continue
            text = " ".join((row.get("protein_description", ""), row.get("gene_names", "")))
            fams = {name for name, pattern in FAMILY_RES if pattern.search(text)}
            for key in (ident, re.sub(r"-mRNA-\d+$", "", ident)):
                families[key] = fams
                descriptions[key] = row.get("protein_description", "")
    return families, descriptions


def main() -> int:
    bpw = load_track_module()
    OUT.mkdir(parents=True, exist_ok=True)
    cands = read_tsv(WORK / "junction_candidates.tsv")

    # unit-pair conservation across PCAM/PPE (all species rows)
    pair_species: dict[str, set[str]] = defaultdict(set)
    for row in cands:
        pair_species[row["source_unit_pair"]].add(row["species"])
    def shared_flag(pair: str) -> str:
        sp = {s.upper() for s in pair_species.get(pair, set())}
        has = ("PCAM" in sp) + ("PPE" in sp)
        return "both" if has == 2 else ("one" if has == 1 else "none")

    window_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    per_label_genes: dict[str, dict[str, list]] = {}
    per_label_families: dict[str, dict[str, set[str]]] = {}
    per_label_desc: dict[str, dict[str, str]] = {}
    track_cache: dict[tuple[str, str], dict[str, object]] = {}

    for species, labels in SPECIES_LABELS.items():
        fam_map, desc_map = functional_families(ROOT, species)
        for label in labels:
            gff = ROOT / "0903" / "karyotype" / "prepared" / label / f"{label}.gff"
            per_label_genes[label] = gene_index(label, gff)
            per_label_families[label] = fam_map
            per_label_desc[label] = desc_map

    for row in cands:
        species = row["species"]
        if species not in SPECIES_LABELS:
            continue
        label = row["target_label"]
        chrom = row["target_chromosome"]
        if label not in per_label_genes or chrom not in per_label_genes[label]:
            continue
        track = "A" if row["subgenome"] == "A" else "B"
        center = int(float(row["candidate_center_bp"]))
        fam_map = per_label_families[label]
        desc_map = per_label_desc[label]
        chrom_len = int(float(row["target_chromosome_length_bp"]))
        gff = ROOT / "0903" / "karyotype" / "prepared" / label / f"{label}.gff"
        genes2 = defaultdict(dict)
        with gff.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[0] != chrom:
                    continue
                gene = re.sub(r"-mRNA-\d+$", "", f[1])
                try:
                    start, end = int(f[2]), int(f[3])
                except ValueError:
                    continue
                old = genes2[gene]
                if not old or (end - start) > (old[1] - old[0]):
                    genes2[gene] = (start, end)
        order = sorted(genes2.items(), key=lambda kv: (kv[1][0] + kv[1][1]) / 2)
        mids = np.array([(v[0] + v[1]) / 2 for _k, v in order])
        fams = [fam_map.get(g, set()) for g, _v in order]
        for width_name, width in WINDOWS.items():
            low, high = max(0, center - width), min(chrom_len, center + width)
            in_win = np.where((mids >= low) & (mids <= high))[0]
            n_genes = int(in_win.size)
            fam_counter: Counter = Counter()
            for idx in in_win:
                for fam in fams[idx]:
                    fam_counter[fam] += 1
                if fams[idx]:
                    gene_rows.append({
                        "species": species, "target_label": label, "track": track,
                        "shared": shared_flag(row["source_unit_pair"]),
                        "candidate_id": row["candidate_id"], "gene": order[idx][0],
                        "families": ";".join(sorted(fams[idx])),
                        "description": desc_map.get(order[idx][0], ""),
                        "window": width_name,
                    })
            te_values: dict[str, float] = {}
            if width_name == "main":
                species_dir = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}[species]
                ident = f"{'A' if track == 'A' else 'B'}{'' if species == 'PPA' else ''}"
                hap = re.sub(r"^p[sz]i?_?", "", label)
                hap_code = label.split("hap")[-1] if "hap" in label else ""
                circos_id = f"{hap_code}_{chrom}"
                key = (species_dir, circos_id)
                if key not in track_cache:
                    cache = {}
                    tdir = ROOT / "0903" / "haploid_repartition" / "circos" / species_dir / "data"
                    for tname in TE_TRACKS:
                        try:
                            parsed = bpw.parse_track(tdir / f"{tname}.txt")
                        except FileNotFoundError:
                            parsed = {}
                        cache[tname] = bpw.TrackIndex(parsed.get(circos_id, []))
                    track_cache[key] = cache
                for tname in TE_TRACKS:
                    te_values[tname] = bpw.TrackIndex.mean(track_cache[key][tname], low, high)
            window_rows.append({
                "species": species, "target_label": label, "track": track,
                "shared": shared_flag(row["source_unit_pair"]),
                "candidate_id": row["candidate_id"], "chromosome": chrom,
                "chrom_len": chrom_len,
                "center_bp": center, "window": width_name, "n_genes": n_genes,
                "n_candidate_genes": int(sum(fam_counter.values())),
                **{f"n_{fam}": fam_counter.get(fam, 0) for fam, _ in CANDIDATES},
                **{f"te_{k}": v for k, v in te_values.items()},
            })

    write_tsv(OUT / "breakpoint_window_metrics.tsv", [r for r in window_rows if r["window"] == "main"])
    write_tsv(OUT / "breakpoint_window_metrics_sensitivity.tsv", window_rows)
    write_tsv(OUT / "breakpoint_window_candidate_genes.tsv", [r for r in gene_rows if r["window"] == "main"])

    # ---------------------------------------------------------------- null
    # chromosome-matched random windows, family-aware gene counts and TE
    main_rows = [r for r in window_rows if r["window"] == "main"]
    by_region: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for r in main_rows:
        by_region[(str(r["species"]), str(r["track"]), str(r["shared"]), str(r["target_label"]) + "|" + str(r["chromosome"]))].append(r)

    # random pools per (label, chrom): gene midpoints + family flags + TE track
    random_pools: dict[tuple[str, str], dict[str, object]] = {}
    for (species, label_chrom) in sorted({(str(r["species"]), str(r["target_label"]) + "|" + str(r["chromosome"])) for r in main_rows}):
        label, chrom = label_chrom.split("|")
        gff = ROOT / "0903" / "karyotype" / "prepared" / label / f"{label}.gff"
        genes = defaultdict(dict)
        with gff.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[0] != chrom:
                    continue
                gene = re.sub(r"-mRNA-\d+$", "", f[1])
                try:
                    start, end = int(f[2]), int(f[3])
                except ValueError:
                    continue
                old = genes[gene]
                if not old or (end - start) > (old[1] - old[0]):
                    genes[gene] = (start, end)
        order = sorted(genes.items(), key=lambda kv: (kv[1][0] + kv[1][1]) / 2)
        mids = np.array([(v[0] + v[1]) / 2 for _k, v in order])
        fam_map = per_label_families[label]
        fam_flags = [1 if fam_map.get(g) else 0 for g, _v in order]
        fam_flags = np.asarray(fam_flags)
        fam_count = np.zeros((len(CANDIDATES), len(order)))
        for fam_i, (fam_name, _p) in enumerate(CANDIDATES):
            fam_count[fam_i] = [1 if fam_name in fam_map.get(g, set()) else 0 for g, _v in order]
        chrom_lens = [int(r["chrom_len"]) for r in main_rows
                      if r["species"] == species and f"{r['target_label']}|{r['chromosome']}" == label_chrom]
        chrom_len = max(chrom_lens) if chrom_lens else (int(mids.max() * 1.05) if len(mids) else 1)
        random_pools[(species, label_chrom)] = {
            "mids": mids, "fam_flags": fam_flags, "fam_count": fam_count,
            "chrom_len": chrom_len, "family_names": [f for f, _ in CANDIDATES],
        }

    enrich_rows: list[dict[str, object]] = []
    family_rows: list[dict[str, object]] = []
    for (species, track, shared, label_chrom), rows in sorted(by_region.items()):
        label, chrom = label_chrom.split("|")
        pool = random_pools.get((species, label_chrom))
        if pool is None or len(rows) == 0:
            continue
        obs_genes = np.array([float(r["n_genes"]) for r in rows])
        obs_cand = np.array([float(r["n_candidate_genes"]) for r in rows])
        reps = N_RANDOM
        null_gene = np.zeros(reps)
        null_cand = np.zeros(reps)
        null_fam = np.zeros((len(CANDIDATES), reps))
        mids = pool["mids"]
        fam_flags = pool["fam_flags"]
        fam_count = pool["fam_count"]
        chrom_len = pool["chrom_len"]
        width = 2 * WINDOWS["main"]  # observed windows are center +/- WINDOWS["main"]
        max_start = max(1, int(chrom_len - width))
        for rep in range(reps):
            starts = RNG.integers(0, max_start, size=len(rows))
            lows = starts
            highs = starts + width
            g_counts = np.zeros(len(rows))
            c_counts = np.zeros(len(rows))
            f_counts = np.zeros((len(CANDIDATES), len(rows)))
            for i in range(len(rows)):
                lo, hi = lows[i], highs[i]
                left = np.searchsorted(mids, lo, side="left")
                right = np.searchsorted(mids, hi, side="right")
                g_counts[i] = right - left
                c_counts[i] = fam_flags[left:right].sum()
                f_counts[:, i] = fam_count[:, left:right].sum(axis=1)
            null_gene[rep] = g_counts.mean()
            null_cand[rep] = c_counts.mean()
            null_fam[:, rep] = f_counts.mean(axis=1)
        obs_gene_mean = obs_genes.mean()
        obs_cand_mean = obs_cand.mean()
        p_gene = (1 + int((null_gene >= obs_gene_mean).sum())) / (reps + 1)
        p_cand = (1 + int((null_cand >= obs_cand_mean).sum())) / (reps + 1)
        enrich_rows.append({
            "species": species, "track": track, "shared": shared,
            "region": label_chrom, "n_windows": len(rows),
            "gene_obs": round(float(obs_gene_mean), 3), "gene_null": round(float(null_gene.mean()), 3),
            "gene_ratio": round(float(obs_gene_mean / null_gene.mean()), 4) if null_gene.mean() else "",
            "gene_p": f"{p_gene:.4g}",
            "candidate_obs": round(float(obs_cand_mean), 3), "candidate_null": round(float(null_cand.mean()), 3),
            "candidate_ratio": round(float(obs_cand_mean / null_cand.mean()), 4) if null_cand.mean() else "",
            "candidate_p": f"{p_cand:.4g}",
        })
        for fam_i, fam_name in enumerate([f for f, _ in CANDIDATES]):
            obs_f = float(np.array([r[f"n_{fam_name}"] for r in rows]).mean())
            null_f = float(null_fam[fam_i].mean())
            p_f = (1 + int((null_fam[fam_i] >= obs_f).sum())) / (reps + 1)
            family_rows.append({
                "species": species, "track": track, "shared": shared,
                "region": label_chrom, "family": fam_name, "n_windows": len(rows),
                "obs_per_window": round(obs_f, 4), "null_per_window": round(null_f, 4),
                "ratio": round(obs_f / null_f, 4) if null_f else "",
                "p": f"{p_f:.4g}",
            })

    # aggregate strata (species pooled) with size-weighted means
    agg_rows: list[dict[str, object]] = []
    strata = defaultdict(list)
    for row in enrich_rows:
        strata[(row["track"], row["shared"])].append(row)
    for (track, shared), rows in sorted(strata.items()):
        n = sum(int(r["n_windows"]) for r in rows)
        gene_obs = sum(float(r["gene_obs"]) * int(r["n_windows"]) for r in rows) / n
        gene_null = sum(float(r["gene_null"]) * int(r["n_windows"]) for r in rows) / n
        cand_obs = sum(float(r["candidate_obs"]) * int(r["n_windows"]) for r in rows) / n
        cand_null = sum(float(r["candidate_null"]) * int(r["n_windows"]) for r in rows) / n
        # combined permutation p: Fisher on region-level p values
        def fisher_combine(ps):
            ps = [max(float(p), 1e-300) for p in ps]
            stat = -2 * sum(math.log(p) for p in ps)
            from scipy.stats import chi2
            return float(chi2.sf(stat, 2 * len(ps)))
        p_gene = fisher_combine([r["gene_p"] for r in rows])
        p_cand = fisher_combine([r["candidate_p"] for r in rows])
        agg_rows.append({
            "track": track, "shared": shared, "n_windows": n, "n_regions": len(rows),
            "gene_obs": round(gene_obs, 3), "gene_null": round(gene_null, 3),
            "gene_ratio": round(gene_obs / gene_null, 4) if gene_null else "",
            "gene_p_combined": f"{p_gene:.4g}",
            "candidate_obs": round(cand_obs, 4), "candidate_null": round(cand_null, 4),
            "candidate_ratio": round(cand_obs / cand_null, 4) if cand_null else "",
            "candidate_p_combined": f"{p_cand:.4g}",
        })
    write_tsv(OUT / "breakpoint_enrichment_by_region.tsv", enrich_rows)
    write_tsv(OUT / "breakpoint_enrichment_by_strata.tsv", agg_rows)
    write_tsv(OUT / "breakpoint_family_by_region.tsv", family_rows)

    # TE stats: per stratum vs chromosome-matched random windows
    te_rows: list[dict[str, object]] = []
    for (species, track, shared, label_chrom), rows in sorted(by_region.items()):
        label, chrom = label_chrom.split("|")
        pool = random_pools.get((species, label_chrom))
        if pool is None:
            continue
        species_dir = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}[species]
        hap_code = label.split("hap")[-1] if "hap" in label else ""
        circos_id = f"{hap_code}_{chrom}"
        tdir = ROOT / "0903" / "haploid_repartition" / "circos" / species_dir / "data"
        for tname in TE_TRACKS:
            try:
                parsed = bpw.parse_track(tdir / f"{tname}.txt")
            except FileNotFoundError:
                continue
            index = bpw.TrackIndex(parsed.get(circos_id, []))
            if index.length <= 0:
                continue
            obs_vals = []
            for r in rows:
                v = index.mean(int(r["center_bp"]) - WINDOWS["main"], int(r["center_bp"]) + WINDOWS["main"])
                if math.isfinite(v):
                    obs_vals.append(v)
            if not obs_vals:
                continue
            obs = float(np.mean(obs_vals))
            null = []
            for _ in range(N_RANDOM):
                start = int(RNG.integers(0, max(1, index.length - 2 * WINDOWS["main"])))
                v = index.mean(start, start + 2 * WINDOWS["main"])
                if math.isfinite(v):
                    null.append(v)
            if not null:
                continue
            null = np.asarray(null)
            p = (1 + int((null >= obs).sum())) / (len(null) + 1)
            te_rows.append({
                "species": species, "track": track, "shared": shared, "region": label_chrom,
                "track_name": tname, "n_windows": len(obs_vals),
                "obs_mean": round(obs, 4), "null_mean": round(float(null.mean()), 4),
                "ratio": round(obs / null.mean(), 4) if null.mean() else "",
                "p": f"{p:.4g}",
            })
    write_tsv(OUT / "breakpoint_te_by_region.tsv", te_rows)

    # strata-level aggregation for TE
    agg_te = defaultdict(list)
    for r in te_rows:
        agg_te[(r["track_name"], r["track"], r["shared"])].append(r)
    te_agg: list[dict[str, object]] = []
    for (tname, track, shared), rows in sorted(agg_te.items()):
        # directional exact binomial on per-region signs (obs>null) plus mean ratio
        ratios = [float(r["ratio"]) for r in rows if r["ratio"] != ""]
        wins = sum(1 for x in ratios if x > 1)
        from scipy.stats import binomtest
        p_sign = binomtest(wins, len(ratios), 0.5).pvalue if ratios else float("nan")
        te_agg.append({
            "track_name": tname, "track": track, "shared": shared, "n_regions": len(rows),
            "mean_ratio": round(float(np.mean(ratios)), 4) if ratios else "",
            "regions_above_null": wins, "sign_p": f"{p_sign:.4g}" if ratios else "",
        })
    write_tsv(OUT / "breakpoint_te_by_strata.tsv", te_agg)

    print("[OK] windows:", len(main_rows), "| regions:", len(enrich_rows), "| TE rows:", len(te_rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
