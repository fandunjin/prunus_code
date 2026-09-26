#!/usr/bin/env python3
"""Ancestry classification of Fig. 5 junction candidates (task 1).

The v4 breakpoint analysis stratified junctions only by whether the same
source-unit signature was seen in PCAM and/or PPE.  This script reclassifies
every racemose junction candidate by cross-species sharing of its WGDI
source/orientation signature cluster (the WGDI-protocol "shared breakpoint =
ancestral event" logic), then repeats the chromosome-matched permutation tests
of the v4 analysis under the new strata:

  * ancestral_diploid_shared : signature recovered in >=2 racemose species and
    in at least one diploid outgroup (PCAM/PPE);
  * ancestral_trio_only     : recovered in >=2 racemose species, no diploid;
  * lineage_specific        : recovered in a single racemose species.

Outputs land in analysis/breakpoint_environment_20260912/ancestry/ and reuse
the frozen window metrics of the v4 run (same candidates, same windows, same
null construction).
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
spec = importlib.util.spec_from_file_location(
    "bpe", ROOT / "scripts" / "breakpoint_environment_analysis.py")
bpe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bpe)

OUT = bpe.OUT / "ancestry"
V4_METRICS = bpe.OUT / "breakpoint_window_metrics.tsv"
TRIO = {"PZI", "PSP", "PPA"}
DIPLOID = {"PCAM", "PPE"}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    bpw = bpe.load_track_module()

    cands = bpe.read_tsv(bpe.WORK / "junction_candidates.tsv")
    clusters = bpe.read_tsv(bpe.WORK / "junction_candidate_clusters.tsv")
    cid2cluster = {}
    for cl in clusters:
        for cid in cl["candidate_ids"].split(","):
            cid2cluster[cid.strip()] = cl

    class_rows = []
    cls_of: dict[str, str] = {}
    for row in cands:
        if row["species"] not in bpe.SPECIES_LABELS:
            continue
        cl = cid2cluster.get(row["candidate_id"], {})
        sp = {s for s in (cl.get("species", "") or "").split(",") if s}
        trio_n = len(sp & TRIO)
        dip_n = len(sp & DIPLOID)
        if trio_n >= 2:
            cls = "ancestral_diploid_shared" if dip_n >= 1 else "ancestral_trio_only"
        else:
            cls = "lineage_specific"
        cls_of[row["candidate_id"]] = cls
        class_rows.append({
            "candidate_id": row["candidate_id"], "species": row["species"],
            "target_label": row["target_label"], "track": "A" if row["subgenome"] == "A" else "B",
            "chromosome": row["target_chromosome"], "source_unit_pair": row["source_unit_pair"],
            "orientation": f"{row['left_orientation']}|{row['right_orientation']}",
            "evidence_class": row["evidence_class"],
            "cluster_id": cl.get("cluster_id", ""), "cluster_species": cl.get("species", ""),
            "trio_species_n": trio_n, "diploid_species_n": dip_n,
            "deep_outgroup": ";".join(sorted(sp & {"GTR", "RCH"})),
            "ancestry_class": cls,
        })
    bpe.write_tsv(OUT / "candidate_ancestry_classes.tsv", class_rows)
    print("classified candidates:", len(class_rows))
    print(Counter(r["ancestry_class"] for r in class_rows))
    print(Counter((r["track"], r["ancestry_class"]) for r in class_rows))

    metrics = bpe.read_tsv(V4_METRICS)
    metrics = [r for r in metrics if r["candidate_id"] in cls_of]
    for r in metrics:
        r["ancestry"] = cls_of[r["candidate_id"]]
    print("main-window rows:", len(metrics))

    # random pools per (species, label|chrom): gene midpoints + family flags
    random_pools: dict[tuple[str, str], dict[str, object]] = {}
    for species, label_chrom in sorted({(r["species"], f"{r['target_label']}|{r['chromosome']}") for r in metrics}):
        label, chrom = label_chrom.split("|")
        gff = bpe.ROOT / "0903" / "karyotype" / "prepared" / label / f"{label}.gff"
        genes: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
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
                old = genes.get(gene)
                if old is None or (end - start) > (old[1] - old[0]):
                    genes[gene] = (start, end)
        order = sorted(genes.items(), key=lambda kv: (kv[1][0] + kv[1][1]) / 2)
        mids = np.array([(v[0] + v[1]) / 2 for _k, v in order])
        fam_map, _desc = bpe.functional_families(bpe.ROOT, species)
        fam_count = np.zeros((len(bpe.CANDIDATES), len(order)))
        for fam_i, (fam_name, _p) in enumerate(bpe.CANDIDATES):
            fam_count[fam_i] = [1 if fam_name in fam_map.get(g, set()) else 0 for g, _v in order]
        chrom_lens = [int(r["chrom_len"]) for r in metrics
                      if r["species"] == species and f"{r['target_label']}|{r['chromosome']}" == label_chrom]
        chrom_len = max(chrom_lens) if chrom_lens else 1
        random_pools[(species, label_chrom)] = {
            "mids": mids, "fam_count": fam_count, "chrom_len": chrom_len,
        }

    by_group: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for r in metrics:
        by_group[(r["species"], r["track"], r["ancestry"],
                  f"{r['target_label']}|{r['chromosome']}")].append(r)

    rg = np.random.default_rng(20260913)
    enrich_rows = []
    family_rows = []
    for (species, track, ancestry, label_chrom), rows in sorted(by_group.items()):
        pool = random_pools.get((species, label_chrom))
        if pool is None or not rows:
            continue
        obs_genes = np.array([float(r["n_genes"]) for r in rows])
        obs_cand = np.array([float(r["n_candidate_genes"]) for r in rows])
        fam_names = [f for f, _ in bpe.CANDIDATES]
        obs_fam = np.array([[float(r[f"n_{f}"]) for r in rows] for f in fam_names])
        reps = bpe.N_RANDOM
        mids, fam_count, chrom_len = pool["mids"], pool["fam_count"], pool["chrom_len"]
        width = 2 * bpe.WINDOWS["main"]
        max_start = max(1, int(chrom_len - width))
        null_gene = np.zeros(reps)
        null_cand = np.zeros(reps)
        null_fam = np.zeros((len(fam_names), reps))
        for rep in range(reps):
            starts = rg.integers(0, max_start, size=len(rows))
            g_counts = np.zeros(len(rows))
            c_counts = np.zeros(len(rows))
            f_counts = np.zeros((len(fam_names), len(rows)))
            for i in range(len(rows)):
                lo, hi = starts[i], starts[i] + width
                left = np.searchsorted(mids, lo, side="left")
                right = np.searchsorted(mids, hi, side="right")
                g_counts[i] = right - left
                c_counts[i] = fam_count[:, left:right].sum(axis=0).clip(0, 1).sum() if right > left else 0.0
                f_counts[:, i] = fam_count[:, left:right].sum(axis=1)
            null_gene[rep] = g_counts.mean()
            null_cand[rep] = c_counts.mean()
            null_fam[:, rep] = f_counts.mean(axis=1)
        p_gene = (1 + int((null_gene >= obs_genes.mean()).sum())) / (reps + 1)
        p_cand = (1 + int((null_cand >= obs_cand.mean()).sum())) / (reps + 1)
        enrich_rows.append({
            "species": species, "track": track, "ancestry": ancestry, "region": label_chrom,
            "n_windows": len(rows),
            "gene_obs": round(float(obs_genes.mean()), 3),
            "gene_null": round(float(null_gene.mean()), 3),
            "gene_ratio": round(float(obs_genes.mean() / null_gene.mean()), 4) if null_gene.mean() else "",
            "gene_p": f"{p_gene:.4g}",
            "candidate_obs": round(float(obs_cand.mean()), 3),
            "candidate_null": round(float(null_cand.mean()), 3),
            "candidate_ratio": round(float(obs_cand.mean() / null_cand.mean()), 4) if null_cand.mean() else "",
            "candidate_p": f"{p_cand:.4g}",
        })
        for fam_i, fam_name in enumerate(fam_names):
            obs_f = float(obs_fam[fam_i].mean())
            null_f = float(null_fam[fam_i].mean())
            p_f = (1 + int((null_fam[fam_i] >= obs_f).sum())) / (reps + 1)
            family_rows.append({
                "species": species, "track": track, "ancestry": ancestry, "region": label_chrom,
                "family": fam_name, "n_windows": len(rows),
                "obs_per_window": round(obs_f, 4), "null_per_window": round(null_f, 4),
                "ratio": round(obs_f / null_f, 4) if null_f else "", "p": f"{p_f:.4g}",
            })

    # aggregate by track x ancestry
    from scipy.stats import chi2, binomtest

    def fisher_combine(ps):
        ps = [max(float(p), 1e-300) for p in ps]
        stat = -2 * sum(math.log(p) for p in ps)
        return float(chi2.sf(stat, 2 * len(ps)))

    agg_rows = []
    strata = defaultdict(list)
    for row in enrich_rows:
        strata[(row["track"], row["ancestry"])].append(row)
    for (track, ancestry), rows in sorted(strata.items()):
        n = sum(int(r["n_windows"]) for r in rows)
        gene_obs = sum(float(r["gene_obs"]) * int(r["n_windows"]) for r in rows) / n
        gene_null = sum(float(r["gene_null"]) * int(r["n_windows"]) for r in rows) / n
        cand_obs = sum(float(r["candidate_obs"]) * int(r["n_windows"]) for r in rows) / n
        cand_null = sum(float(r["candidate_null"]) * int(r["n_windows"]) for r in rows) / n
        agg_rows.append({
            "track": track, "ancestry": ancestry, "n_windows": n, "n_regions": len(rows),
            "gene_obs": round(gene_obs, 3), "gene_null": round(gene_null, 3),
            "gene_ratio": round(gene_obs / gene_null, 4) if gene_null else "",
            "gene_p_combined": f"{fisher_combine([r['gene_p'] for r in rows]):.4g}",
            "candidate_obs": round(cand_obs, 4), "candidate_null": round(cand_null, 4),
            "candidate_ratio": round(cand_obs / cand_null, 4) if cand_null else "",
            "candidate_p_combined": f"{fisher_combine([r['candidate_p'] for r in rows]):.4g}",
        })
    bpe.write_tsv(OUT / "breakpoint_enrichment_by_region.tsv", enrich_rows)
    bpe.write_tsv(OUT / "breakpoint_enrichment_by_ancestry.tsv", agg_rows)

    # family aggregation
    fam_agg = []
    fam_groups = defaultdict(list)
    for row in family_rows:
        fam_groups[(row["track"], row["ancestry"], row["family"])].append(row)
    for (track, ancestry, fam), rows in sorted(fam_groups.items()):
        n = sum(int(r["n_windows"]) for r in rows)
        obs = sum(float(r["obs_per_window"]) * int(r["n_windows"]) for r in rows) / n
        nul = sum(float(r["null_per_window"]) * int(r["n_windows"]) for r in rows) / n
        fam_agg.append({
            "track": track, "ancestry": ancestry, "family": fam, "n_windows": n, "n_regions": len(rows),
            "obs_per_window": round(obs, 4), "null_per_window": round(nul, 4),
            "ratio": round(obs / nul, 4) if nul else "",
            "p_combined": f"{fisher_combine([r['p'] for r in rows]):.4g}",
        })
    bpe.write_tsv(OUT / "breakpoint_family_by_ancestry.tsv", fam_agg)

    # TE by ancestry: reuse frozen per-candidate TE values, rebuild null per pool
    te_track_names = ["te_coverage", "gypsy_coverage", "copia_coverage", "helitron_coverage", "tir_coverage", "gene_density"]
    te_rows = []
    pools_te: dict[tuple[str, str, str], object] = {}
    for (species, track, ancestry, label_chrom), rows in sorted(by_group.items()):
        label, chrom = label_chrom.split("|")
        species_dir = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}[species]
        hap_code = label.split("hap")[-1] if "hap" in label else ""
        circos_id = f"{hap_code}_{chrom}"
        for tname in te_track_names:
            key = (species_dir, circos_id, tname)
            if key not in pools_te:
                tdir = bpe.ROOT / "0903" / "haploid_repartition" / "circos" / species_dir / "data"
                try:
                    parsed = bpw.parse_track(tdir / f"{tname}.txt")
                except FileNotFoundError:
                    pools_te[key] = None
                    continue
                index = bpw.TrackIndex(parsed.get(circos_id, []))
                if index.length <= 0:
                    pools_te[key] = None
                    continue
                null = []
                for _ in range(bpe.N_RANDOM):
                    start = int(rg.integers(0, max(1, index.length - 2 * bpe.WINDOWS["main"])))
                    v = index.mean(start, start + 2 * bpe.WINDOWS["main"])
                    if math.isfinite(v):
                        null.append(v)
                pools_te[key] = (index, np.asarray(null)) if null else None
            entry = pools_te[key]
            if entry is None:
                continue
            index, null = entry
            obs_vals = []
            for r in rows:
                v = index.mean(int(r["center_bp"]) - bpe.WINDOWS["main"], int(r["center_bp"]) + bpe.WINDOWS["main"])
                if math.isfinite(v):
                    obs_vals.append(v)
            if not obs_vals:
                continue
            obs = float(np.mean(obs_vals))
            p = (1 + int((null >= obs).sum())) / (len(null) + 1)
            te_rows.append({
                "species": species, "track": track, "ancestry": ancestry, "region": label_chrom,
                "track_name": tname, "n_windows": len(obs_vals),
                "obs_mean": round(obs, 4), "null_mean": round(float(null.mean()), 4),
                "ratio": round(obs / null.mean(), 4) if null.mean() else "", "p": f"{p:.4g}",
            })
    bpe.write_tsv(OUT / "breakpoint_te_by_region.tsv", te_rows)

    te_agg = []
    te_groups = defaultdict(list)
    for r in te_rows:
        te_groups[(r["track_name"], r["track"], r["ancestry"])].append(r)
    for (tname, track, ancestry), rows in sorted(te_groups.items()):
        ratios = [float(r["ratio"]) for r in rows if r["ratio"] != ""]
        wins = sum(1 for x in ratios if x > 1)
        p_sign = binomtest(wins, len(ratios), 0.5).pvalue if ratios else float("nan")
        te_agg.append({
            "track_name": tname, "track": track, "ancestry": ancestry, "n_regions": len(rows),
            "mean_ratio": round(float(np.mean(ratios)), 4) if ratios else "",
            "regions_above_null": wins,
            "sign_p": f"{p_sign:.4g}" if ratios else "",
        })
    bpe.write_tsv(OUT / "breakpoint_te_by_ancestry.tsv", te_agg)

    print("[OK] groups:", len(enrich_rows), "| family rows:", len(family_rows), "| TE rows:", len(te_rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
