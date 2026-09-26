#!/usr/bin/env python3
"""Pillar-level functional fractionation, homeolog Ka/Ks, and model-based breakpoint analysis (v2).

Uses the fixed A/B physical-track POInT run (#25) per the author's final choice.
Direction 1a  floral-regulator pillar retention bias (A/B) vs genome background
Direction 1b  A-B homeolog Ka/Ks by family; A/B track constraint asymmetry vs outgroups
Direction 3   pillar-state composition around PZI A/B rearrangement junctions
"""
from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, binomtest

ROOT = Path(__file__).resolve().parents[1]
PILLAR = ROOT / "0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping/25_point_fixed_ab_physical_tracks/01_model_run/normalized/point_posterior_summary.tsv"
KAKS_DIR = ROOT / "analysis/wgd_kaks/block_representative_YN"
PAIRS_DIR = KAKS_DIR / "pairs"
KAKS_FILES = KAKS_DIR / "kaks"
BP = ROOT / "0903/karyotype/fig4_breakpoint_candidates"
OUT = ROOT / "analysis/pillar_functional_20260912"
RNG = np.random.default_rng(20260912)
NPERM = 5000

FAMILY_PATTERNS = [
    ("LFY", r"Floricaula|LEAFY protein|\bLEAFY\b"),
    ("AP1", r"Floral homeotic protein APETALA 1|APETALA1\b"),
    ("FUL", r"FRUITFULL|\bFUL\b"),
    ("TFL1/CEN", r"Terminal flower 1\b|TERMINAL FLOWER 1\b"),
    ("FT", r"Flowering locus T|\bFT protein\b"),
    ("SOC1", r"SUPPRESSOR OF OVEREXPRESSION|MADS-box protein SOC1"),
    ("SVP/DAM", r"SHORT VEGETATIVE PHASE|[Dd]ormancy-associated|[Dd]ormancy associated"),
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
FAMILY_RES = [(n, re.compile(p)) for n, p in FAMILY_PATTERNS]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def gk(x: str) -> str:
    x = x.split("__")[-1]
    x = re.sub(r"-mRNA-\d+$", "", x)
    m = re.search(r"(LOC_\d+)", x)
    return m.group(1) if m else x


def build_family_map(species: str) -> dict[str, set[str]]:
    path = ROOT / "0903/annotation/functional" / species / "swissprot_functional_annotation.tsv"
    out: dict[str, set[str]] = defaultdict(set)
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_id = header.index("protein_id")
            i_desc = header.index("protein_description")
            i_gn = header.index("gene_names")
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_id, i_desc, i_gn):
                continue
            text = f[i_desc] + " " + f[i_gn]
            fams = {n for n, r in FAMILY_RES if r.search(text)}
            if fams:
                out[gk(f[i_id])] |= fams
    return out


def bh(pvals: list[float]) -> list[float]:
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    out = [1.0] * len(pvals)
    prev = 1.0
    n = len(pvals)
    for rank, idx in enumerate(reversed(order), 1):
        j = n - rank + 1
        val = min(prev, pvals[idx] * n / j)
        out[idx] = val
        prev = val
    return out


def load_table() -> list[dict[str, str]]:
    return list(csv.DictReader(PILLAR.open(encoding="utf-8"), delimiter="\t"))


A_KEYS = ("PZI_A_gene", "PSP_A_gene", "PPA_A_gene")
B_KEYS = ("PZI_B_gene", "PSP_B_gene", "PPA_B_gene")


def direction_1a(rows):
    fam = {"pzi": build_family_map("pzi"), "psp": build_family_map("psp"), "ppa": build_family_map("ppa")}
    sp_of = {"PZI": "pzi", "PSP": "psp", "PPA": "ppa"}
    pillar_fams = []
    a_ret, b_ret = [], []
    for r in rows:
        fams: set[str] = set()
        for t in A_KEYS + B_KEYS:
            g = r.get(t, "NONE")
            if g and g != "NONE":
                fams |= fam[sp_of[t[:3]]].get(gk(g), set())
        pillar_fams.append(fams)
        a_ret.append(sum(1 for t in A_KEYS if r.get(t, "NONE") != "NONE"))
        b_ret.append(sum(1 for t in B_KEYS if r.get(t, "NONE") != "NONE"))
    a_ret = np.asarray(a_ret, dtype=float)
    b_ret = np.asarray(b_ret, dtype=float)
    print(f"[1a] pillars={len(rows)} mean A={a_ret.mean():.4f} B={b_ret.mean():.4f} (fixed A/B run)")
    rows_out = []
    masks = {"ALL_floral": np.array([bool(f) for f in pillar_fams])}
    for name, _r in FAMILY_RES:
        masks[name] = np.array([name in f for f in pillar_fams])
    bg_delta = b_ret.mean() - a_ret.mean()
    for name, mask in masks.items():
        k = int(mask.sum())
        if k < 8:
            continue
        obs = b_ret[mask].mean() - a_ret[mask].mean()
        idx_pool = np.where((a_ret + b_ret) > 0)[0]
        samples = RNG.choice(idx_pool, size=(NPERM, k), replace=True)
        null = b_ret[samples].mean(axis=1) - a_ret[samples].mean(axis=1)
        p_adv = (1 + int((null >= obs).sum())) / (NPERM + 1)
        p_red = (1 + int((null <= obs).sum())) / (NPERM + 1)
        a_only = int(((a_ret[mask] > 0) & (b_ret[mask] == 0)).sum())
        b_only = int(((b_ret[mask] > 0) & (a_ret[mask] == 0)).sum())
        bg_a = int(((a_ret > 0) & (b_ret == 0)).sum())
        bg_b = int(((b_ret > 0) & (a_ret == 0)).sum())
        orr, pf = fisher_exact([[a_only, b_only], [bg_a, bg_b]])
        rows_out.append({
            "set": name, "n_pillars": k,
            "mean_A_retention": f"{a_ret[mask].mean():.4f}", "mean_B_retention": f"{b_ret[mask].mean():.4f}",
            "B_minus_A": f"{obs:+.4f}", "background_B_minus_A": f"{bg_delta:+.4f}",
            "perm_p_B_advantage": f"{p_adv:.4g}", "perm_p_B_reduction": f"{p_red:.4g}",
            "A_only": a_only, "B_only": b_only, "fisher_p": f"{pf:.4g}",
        })
    ps = [float(x["perm_p_B_advantage"]) for x in rows_out]
    for x, q in zip(rows_out, bh(ps) if ps else []):
        x["q_B_advantage"] = f"{q:.4g}"
    write_tsv(OUT / "pillar_family_retention.tsv", rows_out)
    print("[1a] wrote", len(rows_out), "rows (n>=8 pillars)")


def build_bridge() -> dict[str, str]:
    """canonical LOC -> legacy LOC via exact protein sequence."""
    legacy = {}
    for ph in ("A1", "A2", "B1", "B2"):
        p = ROOT / "0903/haploid_repartition/inputs/pzi/pep" / f"pzi_hap{ph}.pep.fa"
        key = None
        seq = []
        for line in p.open(encoding="utf-8", errors="ignore"):
            if line.startswith(">"):
                if key:
                    legacy.setdefault("".join(seq).replace("*", ""), gk(key))
                key = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if key:
            legacy.setdefault("".join(seq).replace("*", ""), gk(key))
    bridge = {}
    for ph in ("A1", "A2", "B1", "B2"):
        p = ROOT / "analysis/canonical/pzi_hap" / "" / f"pzi_hap{ph}" / "primary.pep.fa"
        p = ROOT / "analysis/canonical" / f"pzi_hap{ph}" / "primary.pep.fa"
        key = None
        seq = []
        for line in p.open(encoding="utf-8", errors="ignore"):
            if line.startswith(">"):
                if key:
                    leg = legacy.get("".join(seq).replace("*", ""))
                    if leg:
                        bridge[gk(key)] = leg
                key = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if key:
            leg = legacy.get("".join(seq).replace("*", ""))
            if leg:
                bridge[gk(key)] = leg
    print(f"[bridge] canonical->legacy exact matches: {len(bridge)}")
    return bridge


def load_kaks(comp: str):
    path = KAKS_FILES / f"{comp}.KaKs_Calculator.YN.tsv"
    out = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["Sequence"]] = row
    return out


def load_pairs(comp: str):
    path = PAIRS_DIR / f"{comp}.alignment_stats.tsv"
    out = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["pair_id"]] = row
    return out


def direction_1b():
    fam = build_family_map("pzi")
    bridge = build_bridge()

    def fams_of(canonical_gene: str) -> set[str]:
        loc = gk(canonical_gene)
        leg = bridge.get(loc, "")
        return fam.get(leg, set()) if leg else set()

    ab_comps = ["pzi_hapA1_vs_pzi_hapB1", "pzi_hapA1_vs_pzi_hapB2",
                "pzi_hapA2_vs_pzi_hapB1", "pzi_hapA2_vs_pzi_hapB2"]
    pair_rows = []
    ka_by_family: dict[str, list[float]] = defaultdict(list)
    ka_bg: list[float] = []
    for comp in ab_comps:
        kaks = load_kaks(comp)
        pairs = load_pairs(comp)
        for pid, p in pairs.items():
            if p.get("status") != "used":
                continue
            k = kaks.get(pid)
            if not k or not k.get("Ka/Ks"):
                continue
            try:
                ratio = float(k["Ka/Ks"])
            except ValueError:
                continue
            fset = fams_of(p.get("query_gene", "")) | fams_of(p.get("ref_gene", ""))
            bridge_ok = (gk(p.get("query_gene", "")) in bridge) or (gk(p.get("ref_gene", "")) in bridge)
            pair_rows.append({"pair_id": pid, "comparison": comp, "Ka/Ks": ratio,
                              "query_gene": gk(p.get("query_gene", "")), "ref_gene": gk(p.get("ref_gene", "")),
                              "bridged": int(bridge_ok), "families": ";".join(sorted(fset))})
            if fset:
                for f in fset:
                    ka_by_family[f].append(ratio)
            else:
                ka_bg.append(ratio)
    write_tsv(OUT / "homeolog_pair_kaks.tsv", pair_rows)
    n_bridge = sum(1 for r in pair_rows if r["bridged"])
    print(f"[1b] A-B used pairs={len(pair_rows)} bridged={n_bridge} bg={len(ka_bg)} floral={sum(len(v) for v in ka_by_family.values())}")
    results = []
    for name in sorted(ka_by_family):
        vals = ka_by_family[name]
        if len(vals) < 5 or len(ka_bg) < 30:
            continue
        u, p = mannwhitneyu(vals, ka_bg, alternative="two-sided")
        results.append({"family": name, "n_pairs": len(vals),
                        "median_KaKs_family": f"{np.median(vals):.4f}",
                        "median_KaKs_background": f"{np.median(ka_bg):.4f}",
                        "ratio": f"{np.median(vals) / np.median(ka_bg):.4f}",
                        "mannwhitney_p": f"{p:.4g}"})
    ps = [float(x["mannwhitney_p"]) for x in results]
    for x, q in zip(results, bh(ps) if ps else []):
        x["q"] = f"{q:.4g}"
    write_tsv(OUT / "homeolog_kaks_by_family.tsv", results)

    # track asymmetry vs outgroups (pair A-left and B-left by shared outgroup gene)
    asym_rows = []
    for out_sp in ("ppe", "pcam"):
        a_pairs = defaultdict(list)
        b_pairs = defaultdict(list)
        for a_comp, store in (("pzi_hapA1", a_pairs), ("pzi_hapB1", b_pairs)):
            comp = f"{a_comp}_vs_{out_sp}"
            kaks = load_kaks(comp)
            pairs = load_pairs(comp)
            for pid, p in pairs.items():
                if p.get("status") != "used":
                    continue
                k = kaks.get(pid)
                if not k:
                    continue
                try:
                    store[p.get("ref_gene", "")].append((float(k["Ka"]), float(k["Ka/Ks"]), gk(p.get("query_gene", ""))))
                except (ValueError, TypeError):
                    continue
        common = [g for g in a_pairs if g in b_pairs and g]
        for g in common:
            ka_a, rat_a, qa = a_pairs[g][0]
            ka_b, rat_b, qb = b_pairs[g][0]
            fset = fams_of(qa) | fams_of(qb)
            asym_rows.append({"outgroup": out_sp, "outgroup_gene": g, "A_gene": qa, "B_gene": qb,
                              "Ka_A": ka_a, "Ka_B": ka_b, "KaKs_A": rat_a, "KaKs_B": rat_b,
                              "diff_A_minus_B": rat_a - rat_b, "families": ";".join(sorted(fset))})
    if asym_rows:
        diffs = np.array([r["diff_A_minus_B"] for r in asym_rows])
        wins_a = int((diffs > 0).sum())
        wins_b = int((diffs < 0).sum())
        p_sign = binomtest(wins_a, wins_a + wins_b, 0.5).pvalue if (wins_a + wins_b) else float("nan")
        write_tsv(OUT / "track_constraint_asymmetry.tsv", asym_rows)
        print(f"[1b-asym] paired loci={len(asym_rows)} A_higher_KaKs={wins_a} B_higher={wins_b} sign p={p_sign:.4g}")
        floral = [r for r in asym_rows if r["families"]]
        if floral:
            print(f"[1b-asym] floral loci={len(floral)} median diff={np.median([r['diff_A_minus_B'] for r in floral]):+.4f}")
    print(f"[1b] family rows={len(results)}; overall A-B median Ka/Ks={np.median([r['Ka/Ks'] for r in pair_rows]):.4f}")


def load_gff_positions(path: Path):
    out = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            g = gk(f[1])
            try:
                mid = (int(f[2]) + int(f[3])) / 2.0
            except ValueError:
                continue
            out.setdefault(g, (f[0], mid))
    return out


def direction_3(rows):
    apos = load_gff_positions(ROOT / "0903/karyotype/prepared/pzi_hapA1/pzi_hapA1.gff")
    bpos = load_gff_positions(ROOT / "0903/karyotype/prepared/pzi_hapB1/pzi_hapB1.gff")
    a_pts, b_pts = [], []
    for r in rows:
        ga, gb = r.get("PZI_A_gene", "NONE"), r.get("PZI_B_gene", "NONE")
        a_present = sum(1 for t in A_KEYS if r.get(t, "NONE") != "NONE")
        b_present = sum(1 for t in B_KEYS if r.get(t, "NONE") != "NONE")
        if a_present + b_present == 0:
            continue
        if ga != "NONE":
            g = gk(ga)
            if g in apos:
                a_pts.append((apos[g][0], apos[g][1], a_present, b_present))
        if gb != "NONE":
            g = gk(gb)
            if g in bpos:
                b_pts.append((bpos[g][0], bpos[g][1], a_present, b_present))
    print(f"[3] mapped pillars A={len(a_pts)} B={len(b_pts)}")
    cands = list(csv.DictReader((BP / "junction_candidates.tsv").open(encoding="utf-8"), delimiter="\t"))
    pair_species = defaultdict(set)
    for c in cands:
        pair_species[c["source_unit_pair"]].add(c["species"].upper())
    half = 200_000
    win_rows = []
    for lane, points in (("pzi_hapA1", a_pts), ("pzi_hapB1", b_pts)):
        for c in cands:
            if c["target_label"] != lane:
                continue
            chrom = c["target_chromosome"]
            center = int(float(c["candidate_center_bp"]))
            st = Counter()
            for cpos, pos, a, b in points:
                if cpos != chrom or abs(pos - center) > half:
                    continue
                st["n"] += 1
                st["A"] += 1 if a > 0 else 0
                st["B"] += 1 if b > 0 else 0
                if a > 0 and b > 0:
                    st["both"] += 1
                elif a > 0:
                    st["A_only"] += 1
                else:
                    st["B_only"] += 1
            if not st["n"]:
                continue
            win_rows.append({"lane": lane, "chromosome": chrom, "candidate_id": c["candidate_id"],
                             "track": "A" if lane.endswith("A1") else "B",
                             "pcam": int("PCAM" in pair_species[c["source_unit_pair"]]),
                             "ppe": int("PPE" in pair_species[c["source_unit_pair"]]),
                             "n_pillars": st["n"], "A_frac": round(st["A"] / st["n"], 4),
                             "B_frac": round(st["B"] / st["n"], 4),
                             "A_only": st["A_only"], "B_only": st["B_only"], "both": st["both"]})
    # stratum summary vs chromosome-matched random windows
    strata_of = {}
    for r in win_rows:
        if r["lane"] == "pzi_hapA1" and r["chromosome"] in ("C1A", "C2A"):
            strata_of[r["candidate_id"]] = "A_path_C1A_C2A"
        elif r["lane"] == "pzi_hapB1" and r["chromosome"] in ("C6B", "C7B"):
            strata_of[r["candidate_id"]] = "B_side_C6B_C7B"
        else:
            strata_of[r["candidate_id"]] = "other"
    out = []
    for s in ("A_path_C1A_C2A", "B_side_C6B_C7B", "other"):
        rs = [r for r in win_rows if strata_of[r["candidate_id"]] == s]
        if not rs:
            continue
        fr = np.array([r["A_only"] + r["B_only"] for r in rs], dtype=float)
        out.append({"stratum": s, "n_windows": len(rs),
                    "A_frac_mean": round(float(np.mean([r["A_frac"] for r in rs])), 4),
                    "B_frac_mean": round(float(np.mean([r["B_frac"] for r in rs])), 4),
                    "A_only_total": int(sum(r["A_only"] for r in rs)),
                    "B_only_total": int(sum(r["B_only"] for r in rs)),
                    "mean_pillars": round(float(np.mean([r["n_pillars"] for r in rs])), 2)})
    write_tsv(OUT / "breakpoint_pillar_windows.tsv", win_rows)
    write_tsv(OUT / "breakpoint_pillar_strata.tsv", out)
    print("[3] strata:", out)


def direction_3b(rows):
    """Non-circular distance analysis: retention of the *other* track near junctions."""
    apos = load_gff_positions(ROOT / "0903/karyotype/prepared/pzi_hapA1/pzi_hapA1.gff")
    bpos = load_gff_positions(ROOT / "0903/karyotype/prepared/pzi_hapB1/pzi_hapB1.gff")
    cands = list(csv.DictReader((BP / "junction_candidates.tsv").open(encoding="utf-8"), delimiter="\t"))
    for lane, posmap, this_set, other_set, chrom_class in (
            ("pzi_hapA1", apos, A_KEYS, B_KEYS, {"C1A": "A_path", "C2A": "A_path"}),
            ("pzi_hapB1", bpos, B_KEYS, A_KEYS, {"C6B": "B_side", "C7B": "B_side"})):
        junctions = defaultdict(list)
        for c in cands:
            if c["target_label"] == lane:
                junctions[c["target_chromosome"]].append(int(float(c["candidate_center_bp"])))
        rows_out = []
        for r in rows:
            g = r.get("PZI_A_gene" if lane == "pzi_hapA1" else "PZI_B_gene", "NONE")
            if g == "NONE":
                continue
            gg = gk(g)
            if gg not in posmap:
                continue
            chrom, mid = posmap[gg]
            if chrom not in junctions:
                continue
            d = min(abs(mid - j) for j in junctions[chrom])
            other_state = sum(1 for t in other_set if r.get(t, "NONE") != "NONE")
            klass = chrom_class.get(chrom, "background")
            rows_out.append({"lane": lane, "chromosome": chrom, "class": klass, "distance_bp": d,
                             "other_track_retention": other_state})
        import numpy as _np
        bins = [(0, 50_000), (50_000, 200_000), (200_000, 1_000_000), (1_000_000, float("inf"))]
        out = []
        for klass in ("A_path", "B_side", "background"):
            sub = [x for x in rows_out if x["class"] == klass]
            if not sub:
                continue
            for lo, hi in bins:
                grp = [x for x in sub if lo <= x["distance_bp"] < hi]
                if not grp:
                    continue
                rate = sum(x["other_track_retention"] for x in grp) / (3 * len(grp))
                out.append({"lane": lane, "class": klass, "distance_bin": f"{lo}-{hi if hi != float('inf') else 'inf'}",
                            "n_pillars": len(grp), "other_track_retention": round(rate, 4)})
        write_tsv(OUT / f"breakpoint_distance_profile_{lane}.tsv", out)
        print(f"[3b] {lane}: profile rows={len(out)}")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_table()
    direction_1a(rows)
    direction_1b()
    direction_3(rows)
    direction_3b(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
