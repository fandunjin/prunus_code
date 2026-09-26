#!/usr/bin/env python3
"""Statistics for the cross-species inflorescence-family analysis.

Writes, under <outdir>/stats/:
  genome_background_genelevel.tsv  gene-level track totals per species
  family_species_counts.tsv        tree-supported copy counts per species/group/track
  family_ab_tests.tsv              per-species Fisher (family vs genome background)
  family_ab_meta.tsv               Mantel-Haenszel pooled OR and Fisher-combined p
  component_fate_tests.tsv         exact McNemar on A-only vs B-only components
  phase_pair_replication.tsv       A1-B1 versus A2-B2 lineage replication (PZI/PSP)
  component_anchor_coverage.tsv    anchored fraction of tree-supported genes
"""
from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import binomtest, chi2, fisher_exact

ROOT = Path(__file__).resolve().parents[1]
PHASE_RE = re.compile(r"^(pzi|psp|ppa|pcam|ppe|gtr)_hap(A1|A2|B1|B2|A|B|S)_")
LOC_RE = re.compile(r"(LOC_\d+)")
TRACK = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "A": "A", "B": "B", "S": "S"}
SPECIES_ORDER = ["pzi", "psp", "ppa"]
GROUP_DISPLAY = [
    ("MADS", "AP1/CAL"), ("MADS", "SEP"), ("MADS", "AG/SHP/STK"), ("MADS", "AGL6"),
    ("MADS", "SOC1"), ("MADS", "SVP/DAM"), ("MADS", "FUL"), ("MADS", "other-MADS"),
    ("PEBP", "FT-like"), ("PEBP", "TFL1-like"), ("PEBP", "MFT-like"), ("PEBP", "PEBP-other"),
    ("C_LFY_FLO", "LFY/FLO"),
    ("SPL", "SPL9/15"), ("SPL", "SPL3/4/5"), ("SPL", "SPL-other"),
    ("TCP", "TCP2/3/4"), ("TCP", "BRC1/2"), ("TCP", "TCP-other"),
    ("KNOX", "STM"), ("KNOX", "KNAT1/2/6"),
]


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genome_background(out: Path) -> dict[tuple[str, str], int]:
    """Gene-level totals per (species, track) from the target protein sets."""
    counts: Counter = Counter()
    seen: set[tuple[str, str, str]] = set()
    with (out / "targets" / "all_species.primary.pep.fa").open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            key = line[1:].split()[0]
            m = PHASE_RE.match(key)
            if not m:
                continue
            loc = LOC_RE.search(key)
            gene = (m.group(1), m.group(2), loc.group(1) if loc else key)
            if gene in seen:
                continue
            seen.add(gene)
            counts[(m.group(1), TRACK[m.group(2)])] += 1
    return dict(counts)


def group_key(row: dict[str, str]) -> tuple[str, str]:
    return row["family"], row["subfamily_group"]


def family_counts(evidence: list[dict[str, str]]) -> Counter:
    counts: Counter = Counter()
    for row in evidence:
        if row["evidence_tier"] != "tree_supported":
            continue
        counts[(row["species"], row["family"], row["subfamily_group"], row["track"])] += 1
    return counts


def run_fisher(counts, bg, out: Path) -> list[dict[str, object]]:
    rows = []
    groups = sorted({(sp, fam, grp) for (sp, fam, grp, _t) in counts})
    for species, family, group in groups:
        a = counts.get((species, family, group, "A"), 0)
        b = counts.get((species, family, group, "B"), 0)
        ta = bg.get((species, "A"), 0)
        tb = bg.get((species, "B"), 0)
        if a + b == 0 or ta == 0 or tb == 0:
            continue
        ca, cb = ta - a, tb - b
        if a == 0 or b == 0 or ca == 0 or cb == 0:
            a2, b2, ca2, cb2 = a + 0.5, b + 0.5, ca + 0.5, cb + 0.5
        else:
            a2, b2, ca2, cb2 = a, b, ca, cb
        odds, p = fisher_exact([[a, b], [ca, cb]], alternative="two-sided")
        log2or = math.log2((a2 * cb2) / (b2 * ca2))
        se = math.sqrt(sum(1.0 / x for x in (a2, b2, ca2, cb2)))
        rows.append({
            "species": species, "family": family, "subfamily_group": group,
            "A_copies": a, "B_copies": b, "A_total": ta, "B_total": tb,
            "log2_odds_ratio": f"{log2or:.4f}", "or": f"{2 ** log2or:.4f}",
            "ci_low": f"{log2or - 1.96 * se:.4f}", "ci_high": f"{log2or + 1.96 * se:.4f}",
            "p_value": f"{p:.6g}", "n": a + b,
        })
    return rows


def run_meta(fisher_rows: list[dict[str, object]], out: Path) -> list[dict[str, object]]:
    by_group: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in fisher_rows:
        by_group[(str(row["family"]), str(row["subfamily_group"]))].append(row)
    rows = []
    for (family, group), items in sorted(by_group.items()):
        pvals = [float(x["p_value"]) for x in items]
        k = len(pvals)
        stat = -2.0 * sum(math.log(max(p, 1e-300)) for p in pvals)
        combined_p = float(chi2.sf(stat, 2 * k))
        R = S = P = Q = 0.0
        for x in items:
            a = float(x["A_copies"])
            b = float(x["B_copies"])
            c = float(x["A_total"]) - a
            d = float(x["B_total"]) - b
            if min(a, b, c, d) == 0:
                a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            n = a + b + c + d
            r = a * d / n
            s = b * c / n
            R += r
            S += s
            P += (a + d) / n * r
            Q += (b + c) / n * s
        if R > 0 and S > 0:
            mh_or = R / S
            var = (P + Q) / (2 * R * S)
        else:
            mh_or = float("nan")
            var = float("nan")
        log2_mh = math.log2(mh_or) if mh_or == mh_or and mh_or > 0 else float("nan")
        se = math.sqrt(var) / math.log(2) if var == var else float("nan")
        directions = [1 if float(x["log2_odds_ratio"]) > 0 else -1 for x in items]
        n_pos = sum(1 for x in directions if x > 0)
        rows.append({
            "family": family, "subfamily_group": group, "species_n": k,
            "species_tested": ";".join(str(x["species"]) for x in items),
            "mh_or": f"{mh_or:.4f}", "mh_log2_or": f"{log2_mh:.4f}",
            "mh_ci_low": f"{log2_mh - 1.96 * se:.4f}", "mh_ci_high": f"{log2_mh + 1.96 * se:.4f}",
            "fisher_combined_p": f"{combined_p:.6g}",
            "direction_A_species": n_pos,
            "direction_consistency": f"{max(n_pos, k - n_pos)}/{k}",
        })
    return rows


def run_components(out: Path, components: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fate_rows: list[dict[str, object]] = []
    by_species_family: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for comp in components:
        species = comp["species"]
        fate = comp["copy_fate"]
        groups = [g for g in comp["subfamily_groups"].split(";") if g]
        for group in groups:
            # attribute component to the subfamily group of its supported members
            by_species_family[(species, group)][fate] += 1
    for (species, group), counter in sorted(by_species_family.items()):
        aonly = counter.get("A_only", 0)
        bonly = counter.get("B_only", 0)
        both = counter.get("A_and_B", 0)
        n_inf = aonly + bonly
        if n_inf > 0:
            p = float(binomtest(aonly, n_inf, 0.5).pvalue)
        else:
            p = float("nan")
        fate_rows.append({
            "species": species, "subfamily_group": group,
            "components_A_only": aonly, "components_B_only": bonly, "components_A_and_B": both,
            "n_informative": n_inf, "p_exact_sign": f"{p:.6g}" if n_inf else "",
            "direction": "A>B" if aonly > bonly else "B>A" if bonly > aonly else "tie",
        })
    pair_rows: list[dict[str, object]] = []
    by_species_pair: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for comp in components:
        species = comp["species"]
        pairs = [p for p in comp["anchor_pairs"].split(";") if p]
        for pair in pairs:
            for group in [g for g in comp["subfamily_groups"].split(";") if g]:
                by_species_pair[(species, f"{pair}|{group}")][comp["copy_fate"]] += 1
    for key, counter in sorted(by_species_pair.items()):
        species, rest = key
        pair, group = rest.split("|", 1)
        pair_rows.append({
            "species": species, "phase_pair": pair, "subfamily_group": group,
            "A_only": counter.get("A_only", 0), "B_only": counter.get("B_only", 0),
            "A_and_B": counter.get("A_and_B", 0),
        })
    return fate_rows, pair_rows


def gene_key_of(record_id: str) -> str:
    return re.sub(r"-mRNA-\d+$", "", record_id)


def run_coverage(evidence: list[dict[str, str]], edges: list[dict[str, str]]) -> list[dict[str, object]]:
    anchored: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        species = e["species"]
        if e["left_in_target"]:
            anchored[species].add(gene_key_of(str(e["left_id"])))
        if e["right_in_target"]:
            anchored[species].add(gene_key_of(str(e["right_id"])))
    counts: Counter = Counter()
    anchored_counts: Counter = Counter()
    for row in evidence:
        if row["evidence_tier"] != "tree_supported":
            continue
        key = (row["species"], row["family"], row["subfamily_group"])
        counts[key] += 1
        if row["gene_key"] in anchored.get(row["species"], set()):
            anchored_counts[key] += 1
    rows = []
    for key, n in sorted(counts.items()):
        species, family, group = key
        a = anchored_counts[key]
        rows.append({
            "species": species, "family": family, "subfamily_group": group,
            "tree_supported": n, "anchored": a, "anchor_rate": f"{a / n:.4f}",
        })
    return rows


def run_stats(out: Path, evidence: list[dict[str, str]], edges: list[dict[str, str]], logfile: Path) -> None:
    stats_dir = out / "stats"
    bg = genome_background(out)
    write_rows(stats_dir / "genome_background_genelevel.tsv",
               [{"species": s, "track": t, "n": n} for (s, t), n in sorted(bg.items())])
    counts = family_counts(evidence)
    write_rows(stats_dir / "family_species_counts.tsv",
               [{"species": s, "family": f, "subfamily_group": g, "track": tr, "n": n}
                for (s, f, g, tr), n in sorted(counts.items())])
    fisher_rows = run_fisher(counts, bg, out)
    write_rows(stats_dir / "family_ab_tests.tsv", fisher_rows)
    meta_rows = run_meta(fisher_rows, out)
    write_rows(stats_dir / "family_ab_meta.tsv", meta_rows)
    comp_path = out / "copy_fate" / "components_crossspecies.tsv"
    components = read_table(comp_path) if comp_path.exists() else []
    fate_rows, pair_rows = run_components(out, components)
    write_rows(stats_dir / "component_fate_tests.tsv", fate_rows)
    write_rows(stats_dir / "phase_pair_replication.tsv", pair_rows)
    coverage_rows = run_coverage(evidence, edges)
    write_rows(stats_dir / "component_anchor_coverage.tsv", coverage_rows)
    with logfile.open("a", encoding="utf-8") as fh:
        fh.write(f"[stats] fisher rows={len(fisher_rows)} meta rows={len(meta_rows)} fate rows={len(fate_rows)}\n")
    print(f"[stats] fisher={len(fisher_rows)} meta={len(meta_rows)} fate={len(fate_rows)}", flush=True)
