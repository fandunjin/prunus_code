#!/usr/bin/env python3
"""KNOX/BELL breakpoint-proximity review for Fig. 5 (task 2).

Inputs are the frozen v4 outputs of breakpoint_environment_analysis.py plus the
WGDI junction/curation tables. For every KNOX/BELL hit candidate this review:

  1. collapses hits into their source/orientation transition clusters;
  2. measures the KNOX/BELL local neighbourhood (200 kb / 500 kb / 1 Mb) and
     the genome-wide KNOX/BELL cluster (tandem-array hypothesis);
  3. extracts the local WGDI collinearity for the PZI phase sets around each
     junction, to show whether the flanks map to distinct reference intervals
     or to a single unit with an orientation flip (inversion vs boundary
     attraction enriches the raw plotting table only);
  4. writes locus-level tables for the review figure.

Outputs land in analysis/knox_breakpoint_review_20260912/.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "0903" / "karyotype" / "fig4_breakpoint_candidates"
PREP = ROOT / "0903" / "karyotype" / "prepared"
FUNC = ROOT / "0903" / "annotation" / "functional"
V4 = ROOT / "analysis" / "breakpoint_environment_20260912"
OUT = ROOT / "analysis" / "knox_breakpoint_review_20260912"

KNOX_RE = re.compile(r"knotted-1-like|Homeobox protein KNAT|KNOTTED|BELL homeodomain", re.I)
RADII = {"200kb": 200_000, "500kb": 500_000, "1Mb": 1_000_000}
ARRAY_GAP = 100_000  # tandem array = neighbouring KNOX/BELL genes <100 kb apart


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def load_genes(label: str) -> dict[str, list[tuple[float, int, int, str]]]:
    """chromosome -> sorted gene rows (mid, start, end, gene)."""
    gff = PREP / label / f"{label}.gff"
    genes: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    with gff.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            gene = re.sub(r"-mRNA-\d+$", "", f[1])
            try:
                start, end = int(f[2]), int(f[3])
            except ValueError:
                continue
            old = genes[f[0]].get(gene)
            if old is None or (end - start) > (old[1] - old[0]):
                genes[f[0]][gene] = (start, end)
    out: dict[str, list[tuple[float, int, int, str]]] = {}
    for chrom, gmap in genes.items():
        rows = [((s + e) / 2.0, s, e, g) for g, (s, e) in gmap.items()]
        rows.sort()
        out[chrom] = rows
    return out


def load_knox_species(species: str) -> dict[str, str]:
    """gene -> description for KNOX/BELL genes of one species."""
    path = FUNC / species.lower() / "swissprot_functional_annotation.tsv"
    hits: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ident = row.get("protein_id", "")
            if not ident:
                continue
            text = " ".join((row.get("protein_description", ""), row.get("gene_names", "")))
            if KNOX_RE.search(text):
                key = re.sub(r"-mRNA-\d+$", "", ident)
                hits[key] = row.get("protein_description", "")
    return hits


def genomewide_clusters(knox: dict[str, str], genes: dict[str, list[tuple[float, int, int, str]]]):
    rows = []
    for chrom, gmap in genes.items():
        mid_by_gene = {g: mid for mid, _s, _e, g in gmap}
        knox_mid = sorted((mid_by_gene[g], g) for g in knox if g in mid_by_gene)
        if not knox_mid:
            continue
        cluster = [knox_mid[0]]
        for mid, g in knox_mid[1:]:
            if mid - cluster[-1][0] <= ARRAY_GAP:
                cluster.append((mid, g))
            else:
                rows.append((chrom, cluster))
                cluster = [(mid, g)]
        rows.append((chrom, cluster))
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- KNOX hit candidates from the frozen v4 output
    v4_genes = read_tsv(V4 / "breakpoint_window_candidate_genes.tsv")
    knox_hit_ids = sorted({r["candidate_id"] for r in v4_genes if "KNOX/BELL" in r["families"]})
    print("KNOX/BELL hit candidates:", len(knox_hit_ids))

    cands = read_tsv(WORK / "junction_candidates.tsv")
    cand_by_id = {r["candidate_id"]: r for r in cands}
    v4_hits = [r for r in v4_genes if r["candidate_id"] in knox_hit_ids]

    clusters = read_tsv(WORK / "junction_candidate_clusters.tsv")
    cid2cluster = {}
    for cl in clusters:
        for cid in cl["candidate_ids"].split(","):
            cid2cluster[cid.strip()] = cl

    # ---- per species/label caches
    gene_cache: dict[str, dict[str, list[tuple[float, int, int, str]]]] = {}
    knox_cache: dict[str, dict[str, str]] = {}
    species_labels: dict[str, set[str]] = defaultdict(set)
    for label, species in sorted({(cand_by_id[cid]["target_label"], cand_by_id[cid]["species"]) for cid in knox_hit_ids}):
        gene_cache[label] = load_genes(label)
        species_labels[species].add(label)
        if species.lower() not in knox_cache and (FUNC / species.lower()).is_dir():
            knox_cache[species.lower()] = load_knox_species(species.lower())

    # ---- locus neighbourhoods
    locus_rows = []
    for row in v4_hits:
        cand = cand_by_id[row["candidate_id"]]
        label, species, chrom = cand["target_label"], cand["species"], cand["target_chromosome"]
        center = int(float(cand["candidate_center_bp"]))
        genes = gene_cache.get(label, {})
        knox = knox_cache.get(species.lower(), {})
        all_genes = genes.get(chrom, [])
        for radius_name, radius in RADII.items():
            in_win = [r for r in all_genes if abs(r[0] - center) <= radius]
            knox_in = [(abs(r[0] - center), r[3]) for r in in_win if r[3] in knox]
            nearest = min((d for d, _g in knox_in), default=None)
            locus_rows.append({
                "candidate_id": cand["candidate_id"], "species": species, "target_label": label,
                "track": "A" if cand["subgenome"] == "A" else "B", "chromosome": chrom,
                "junction_center_bp": center, "radius": radius_name,
                "n_genes_window": len(in_win),
                "n_knox_window": len(knox_in),
                "knox_genes": ";".join(sorted(g for _d, g in knox_in)),
                "nearest_knox_distance_bp": nearest if nearest is not None else "",
                "source_unit_pair": cand["source_unit_pair"],
                "orientation": f"{cand['left_orientation']}|{cand['right_orientation']}",
                "evidence_class": cand["evidence_class"],
                "bracket_width_bp": cand["candidate_bracket_width_bp"],
                "interblock_gap_genes": cand["interblock_gap_genes"],
                "cluster_id": cid2cluster.get(cand["candidate_id"], {}).get("cluster_id", ""),
                "cluster_species": cid2cluster.get(cand["candidate_id"], {}).get("species", ""),
                "cluster_targets_n": cid2cluster.get(cand["candidate_id"], {}).get("target_count", ""),
            })
    write_tsv(OUT / "knox_locus_neighbourhood.tsv", locus_rows)

    # ---- genome-wide KNOX clusters
    cluster_rows = []
    for species in sorted(knox_cache):
        knox = knox_cache[species]
        for label in sorted(species_labels.get(species.upper(), set())):
            genes = gene_cache[label]
            for ci, (chrom, cluster) in enumerate(genomewide_clusters(knox, genes), start=1):
                mids = [mid for mid, _g in cluster]
                cluster_rows.append({
                    "species": species, "target_label": label, "chromosome": chrom,
                    "cluster_id": f"{label}_{chrom}_K{ci}", "n_knox_genes": len(cluster),
                    "start_bp": int(min(mids)), "end_bp": int(max(mids)),
                    "span_bp": int(max(mids) - min(mids)),
                    "genes": ";".join(g for _m, g in cluster),
                })
    write_tsv(OUT / "knox_genomewide_clusters.tsv", cluster_rows)

    # ---- local WGDI collinearity around the PZI loci
    coll_rows = []
    for label in ("pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2"):
        coll = WORK / "mapping" / label / f"prunus9_vs_{label}.collinearity"
        if not coll.is_file():
            continue
        genes = gene_cache.get(label) or load_genes(label)
        gene_cache[label] = genes
        pos = {}
        for chrom, gmap in genes.items():
            for mid, _s, _e, g in gmap:
                pos[g] = (chrom, mid)
        cands_label = [r for r in v4_hits if cand_by_id[r["candidate_id"]]["target_label"] == label]
        for row in cands_label:
            cand = cand_by_id[row["candidate_id"]]
            center = int(float(cand["candidate_center_bp"]))
            chrom = cand["target_chromosome"]
            unit = {}
            with coll.open(encoding="utf-8", errors="replace") as fh:
                header_units = ""
                for line in fh:
                    if line.startswith("#"):
                        header_units = line.strip().lstrip("#").strip()
                        continue
                    f = line.split()
                    if len(f) < 5:
                        continue
                    ref_gene, ref_order, tgt_gene, tgt_order, orient = f[0], f[1], f[2], f[3], f[4]
                    loc = pos.get(re.sub(r"-mRNA-\d+$", "", tgt_gene).split("__")[-1]) or pos.get(tgt_gene)
                    if loc is None or loc[0] != chrom:
                        continue
                    mid = loc[1]
                    if abs(mid - center) > 2_000_000:
                        continue
                    ref_unit = ref_gene.split("__")[0] if "__" in ref_gene else "?"
                    coll_rows.append({
                        "candidate_id": cand["candidate_id"], "target_label": label,
                        "track": "A" if cand["subgenome"] == "A" else "B", "chromosome": chrom,
                        "junction_center_bp": center, "target_gene_bp": int(mid),
                        "target_gene": tgt_gene, "ref_unit": ref_unit, "ref_order": ref_order,
                        "target_order": tgt_order, "orientation": orient, "alignment": header_units,
                        "knox": "yes" if re.sub(r"-mRNA-\d+$", "", tgt_gene) in knox_cache.get("pzi", {}) else "no",
                    })
    write_tsv(OUT / "knox_collinearity_local.tsv", coll_rows)

    print("locus rows:", len(locus_rows), "| genomewide clusters:", len(cluster_rows), "| collinearity rows:", len(coll_rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
