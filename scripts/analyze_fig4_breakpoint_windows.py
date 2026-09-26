#!/usr/bin/env python3
"""Analyse 0903 WGDI junction candidates around genes and repeat tracks.

The input candidate table is deliberately neutral: only source-unit and
orientation transitions are accepted.  Event names are taken from the
curation table only when the user has independently assigned them.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SPECIES = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
TRACKS = ("gene_density", "te_coverage", "gypsy_coverage")
TRACK_FILES = {track: f"{track}.txt" for track in TRACKS}


class TrackIndex:
    """Prefix-integral index for fast interval means on a CircOS track."""

    def __init__(self, rows: list[tuple[int, int, float]]) -> None:
        rows = sorted(rows)
        self.starts = np.asarray([row[0] for row in rows], dtype=np.int64)
        self.ends = np.asarray([row[1] for row in rows], dtype=np.int64)
        self.values = np.asarray([row[2] for row in rows], dtype=float)
        widths = self.ends - self.starts
        self.integral = np.concatenate(([0.0], np.cumsum(widths * self.values)))
        self.coverage = np.concatenate(([0], np.cumsum(widths)))
        self.length = int(self.ends[-1]) if len(self.ends) else 0

    def mean(self, low: int, high: int) -> float:
        if high <= low or not len(self.starts):
            return float("nan")
        left = int(np.searchsorted(self.ends, low, side="right"))
        right = int(np.searchsorted(self.starts, high, side="left")) - 1
        if left > right or left >= len(self.starts) or right < 0:
            return float("nan")
        left = max(left, 0)
        right = min(right, len(self.starts) - 1)
        weighted = float(self.integral[right + 1] - self.integral[left])
        covered = int(self.coverage[right + 1] - self.coverage[left])
        if self.starts[left] < low:
            delta = min(low, int(self.ends[left])) - int(self.starts[left])
            if delta > 0:
                weighted -= delta * float(self.values[left])
                covered -= delta
        if self.ends[right] > high:
            delta = int(self.ends[right]) - max(high, int(self.starts[right]))
            if delta > 0:
                weighted -= delta * float(self.values[right])
                covered -= delta
        return weighted / covered if covered > 0 else float("nan")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "svg.fonttype": "none", "axes.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_track(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    result: dict[str, list[tuple[int, int, float]]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) < 4:
                continue
            try:
                start, end, value = int(float(fields[1])), int(float(fields[2])), float(fields[3])
            except ValueError:
                continue
            result.setdefault(fields[0], []).append((start, end, value))
    for key in result:
        result[key].sort()
    return result


def target_track_id(label: str, chromosome: str) -> tuple[str, str] | None:
    match = re.match(r"^(pzi|psp)_(hapA1|hapA2|hapB1|hapB2)$", label)
    if match:
        hap = match.group(2).replace("hap", "")
        return match.group(1), f"{hap}_{chromosome}"
    match = re.match(r"^(ppa)_(hapA|hapB)$", label)
    if match:
        hap = match.group(2).replace("hap", "")
        return match.group(1), f"{hap}_{chromosome}"
    return None


def bh_fdr(pvalues: list[float]) -> list[float]:
    indexed = sorted(enumerate(pvalues), key=lambda item: item[1])
    output = [1.0] * len(pvalues)
    previous = 1.0
    n = len(pvalues)
    for rank, (index, value) in reversed(list(enumerate(indexed, 1))):
        adjusted = min(previous, value * n / rank)
        output[index] = adjusted
        previous = adjusted
    return output


def bootstrap_ratio(case: list[float], control: list[float], rng: random.Random, rounds: int = 1000) -> tuple[float, float, float]:
    if not case or not control:
        return float("nan"), float("nan"), float("nan")
    case_mean = float(np.mean(case))
    control_mean = float(np.mean(control))
    ratio = case_mean / control_mean if control_mean != 0 else float("inf")
    values: list[float] = []
    for _ in range(rounds):
        a = [case[rng.randrange(len(case))] for _ in case]
        b = [control[rng.randrange(len(control))] for _ in control]
        denominator = float(np.mean(b))
        values.append(float(np.mean(a)) / denominator if denominator else float("inf"))
    low, high = np.percentile(values, [2.5, 97.5])
    return ratio, float(low), float(high)


def random_window_values(track: TrackIndex, length: int, count: int, rng: random.Random) -> list[float]:
    if track.length <= length:
        return []
    max_start = track.length - length
    values = []
    for _ in range(count):
        start = rng.randint(0, max_start)
        values.append(track.mean(start, start + length))
    return [value for value in values if math.isfinite(value)]


def functional_map(root: Path, species: str) -> dict[str, dict[str, str]]:
    path = root / "annotation" / "functional" / species / "swissprot_functional_annotation.tsv"
    if not path.is_file():
        return {}
    result: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ident = row.get("protein_id", "") or row.get("qseqid", "")
            if ident:
                result[ident] = row
                result[ident.removesuffix("-mRNA-1")] = row
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="0903 directory")
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    parser.add_argument("--random-windows", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260908)
    args = parser.parse_args()
    if args.window_bp <= 0 or args.random_windows < 100:
        raise ValueError("window-bp must be positive and random-windows must be >=100")

    root = args.root.resolve()
    work = root / "karyotype" / "fig4_breakpoint_candidates"
    candidates = read_tsv(work / "junction_candidates.tsv")
    rng = random.Random(args.seed)
    track_cache: dict[tuple[str, str], dict[str, TrackIndex]] = {}
    c_rows: list[dict[str, object]] = []
    density_values: dict[tuple[str, str], list[float]] = {}

    for row in candidates:
        if row.get("species") not in SPECIES:
            continue
        mapped = target_track_id(row["target_label"], row["target_chromosome"])
        if mapped is None:
            continue
        species_dir, ident = mapped
        track_dir = root / "haploid_repartition" / "circos" / species_dir / "data"
        cache_key = (species_dir, ident)
        if cache_key not in track_cache:
            track_cache[cache_key] = {}
            for track in TRACKS:
                parsed = parse_track(track_dir / TRACK_FILES[track])
                # Each target chromosome has one CircOS track; retain an empty
                # index when a file lacks the mapped chromosome identifier.
                track_cache[cache_key][track] = TrackIndex(parsed.get(ident, []))
        center = int(float(row["candidate_center_bp"]))
        low = max(0, center - args.window_bp // 2)
        high = center + args.window_bp // 2
        for track in TRACKS:
            value = track_cache[cache_key][track].mean(low, high)
            if math.isfinite(value):
                density_values.setdefault((species_dir, track), []).append(value)
                c_rows.append({
                    "candidate_id": row["candidate_id"], "target_label": row["target_label"], "species": row["species"],
                    "subgenome": row["subgenome"], "target_chromosome": row["target_chromosome"], "circos_id": ident,
                    "evidence_class": row["evidence_class"], "source_unit_pair": row["source_unit_pair"], "window_low_bp": low,
                    "window_high_bp": high, "track": track, "density": value,
                })

    c_fields = ["candidate_id", "target_label", "species", "subgenome", "target_chromosome", "circos_id", "evidence_class", "source_unit_pair", "window_low_bp", "window_high_bp", "track", "density"]
    with (work / "breakpoint_window_densities.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=c_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(c_rows)

    # Enrichment compares each candidate window with same-track random windows
    # on the same chromosome. This is a density ratio, not a count enrichment.
    by_key: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    for row in c_rows:
        by_key.setdefault((str(row["target_label"]), str(row["circos_id"]), str(row["track"])), []).append(row)
    enrich_rows: list[dict[str, object]] = []
    for (label, ident, track), rows in by_key.items():
        species_dir = target_track_id(label, str(rows[0]["target_chromosome"]))[0]
        source_track = track_cache[(species_dir, ident)][track]
        case = [float(row["density"]) for row in rows]
        controls: list[float] = []
        control_means: list[float] = []
        for row in rows:
            length = int(row["window_high_bp"]) - int(row["window_low_bp"])
            samples = random_window_values(source_track, length, args.random_windows, rng)
            controls.extend(samples)
            if samples:
                control_means.append(float(np.mean(samples)))
        ratio, low, high = bootstrap_ratio(case, control_means, rng)
        observed = float(np.mean(case))
        null_means = np.asarray(controls)
        if len(null_means):
            pvalue = (1 + int(np.sum(null_means >= observed))) / (len(null_means) + 1)
        else:
            pvalue = float("nan")
        enrich_rows.append({
            "target_label": label, "circos_id": ident, "track": track, "n_breakpoint_windows": len(case),
            "n_random_windows": len(controls), "breakpoint_mean_density": observed,
            "random_mean_density": float(np.mean(controls)) if controls else float("nan"), "enrichment_ratio": ratio,
            "ci95_low": low, "ci95_high": high, "permutation_p": pvalue,
        })
    pvalues = [float(row["permutation_p"]) if math.isfinite(float(row["permutation_p"])) else 1.0 for row in enrich_rows]
    fdr = bh_fdr(pvalues)
    for row, value in zip(enrich_rows, fdr):
        row["fdr_bh"] = value
    enrich_fields = ["target_label", "circos_id", "track", "n_breakpoint_windows", "n_random_windows", "breakpoint_mean_density", "random_mean_density", "enrichment_ratio", "ci95_low", "ci95_high", "permutation_p", "fdr_bh"]
    with (work / "breakpoint_window_enrichment.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=enrich_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(enrich_rows)

    # Candidate genes and lightweight floral-development screening.
    keywords = re.compile(r"flower|floral|inflorescence|meristem|APETALA|SOC1|LFY|TFL1|SEP|AGAMOUS|AP1|FLC|FT|SVP|WUS|CLV|TCP", re.I)
    gene_rows: list[dict[str, object]] = []
    for row in candidates:
        species_code = SPECIES.get(row.get("species", ""), "")
        if not species_code:
            continue
        annotations = functional_map(root, species_code)
        for side in ("left", "right"):
            gene = row.get(f"{side}_anchor_gene", "")
            annotation = annotations.get(gene, annotations.get(gene.removesuffix("-mRNA-1"), {}))
            description = annotation.get("protein_description", "") if annotation else ""
            gene_names = annotation.get("gene_names", "") if annotation else ""
            text = f"{description} {gene_names}"
            gene_rows.append({
                "candidate_id": row["candidate_id"], "target_label": row["target_label"], "species": row["species"],
                "subgenome": row["subgenome"], "target_chromosome": row["target_chromosome"], "side": side,
                "gene_id": gene, "protein_description": description, "gene_names": gene_names,
                "floral_development_keyword_match": bool(keywords.search(text)),
                "functional_annotation_status": "matched" if annotation else "not_found",
                "orthogroup_status": "requires_LOC_to_orthogroup_ID_mapping",
            })
    gene_fields = ["candidate_id", "target_label", "species", "subgenome", "target_chromosome", "side", "gene_id", "protein_description", "gene_names", "floral_development_keyword_match", "functional_annotation_status", "orthogroup_status"]
    with (work / "breakpoint_anchor_genes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(gene_rows)

    (work / "analysis_methods.txt").write_text(
        f"Breakpoint window width: {args.window_bp} bp centered on the candidate interval midpoint.\n"
        "Candidate windows are derived from WGDI source-unit/orientation transitions and nearest flanking anchor genes.\n"
        "Random windows are sampled on the same CircOS chromosome and use the same width.\n"
        "Enrichment ratio = mean breakpoint-window density / mean matched-random-window density.\n"
        "95% CI = 1,000-round non-parametric bootstrap over candidate and matched-random window means; permutation P compares each observed density to the random-window null.\n"
        "FDR = Benjamini-Hochberg across all target-label by track tests.\n"
        "Event labels RCT, EEJ and inversion are not inferred automatically; they must be assigned from the curation template after reciprocal and outgroup review.\n"
        "Candidate anchor genes are screened against the 0903 Swiss-Prot functional tables. LOC-to-Orthogroup mapping remains required for gene-family claims.\n",
        encoding="utf-8",
    )
    print(f"[OK] densities: {len(c_rows):,}")
    print(f"[OK] enrichment tests: {len(enrich_rows):,}")
    print(f"[OK] anchor genes: {len(gene_rows):,}")


if __name__ == "__main__":
    main()
