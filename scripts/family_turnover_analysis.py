#!/usr/bin/env python3
"""Gene-family turnover on the racemose lineage (OrthoFinder OGs x species tree).

- Species-level copy numbers from Orthogroups.txt (trio = A+B).
- Simple post-WGD retention ratio: trio copies / (2 x outgroup median).
- Tree-based squared-change (least-squares) ancestral states on a curated
  seven-lineage tree: (rch,(gtr,((pcam,ppe),(ppa,(psp,pzi))))).
- Floral annotation of OGs through exact-protein bridging of the six-species
  candidate records to the frozen OrthoFinder gene IDs.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
OF = ROOT / "analysisphylogenetics/02_orthofinder/main/Results_Aug25/Orthogroups/Orthogroups.txt"
FROZEN = ROOT / "analysisphylogenetics/01_frozen_inputs/main/pep"
CAND = ROOT / "analysis/inflorescence_family_pipeline_crossspecies6_20260912/candidates_crossspecies.tsv"
OUT = ROOT / "analysis/family_turnover_20260912"

TRACKS = ["pzi_A", "pzi_B", "psp_A", "psp_B", "ppa_A", "ppa_B", "pcam", "ppe", "gtr", "rch",
          "pdu", "pav", "psal", "pmal", "pmu", "pspe", "par"]
SPECIES = {"pzi": ["pzi_A", "pzi_B"], "psp": ["psp_A", "psp_B"], "ppa": ["ppa_A", "ppa_B"],
           "pcam": ["pcam"], "ppe": ["ppe"], "gtr": ["gtr"], "rch": ["rch"]}
TREE = "(rch,(gtr,((pcam,ppe),(ppa,(psp,pzi)))));"
CAND_PEP = {
    "pzi": [ROOT / "analysis/canonical" / f"pzi_hap{p}" / "primary.pep.fa" for p in ("A1", "A2", "B1", "B2")],
    "psp": [ROOT / "0903/haploid_repartition/inputs/psp/pep" / f"psp_hap{p}.pep.fa" for p in ("A1", "A2", "B1", "B2")],
    "ppa": [ROOT / "0903/haploid_repartition/inputs/ppa/pep" / f"ppa_hap{p}.pep.fa" for p in ("A", "B")],
    "pcam": [ROOT / "analysis/canonical/pcam/primary.pep.fa"],
    "ppe": [ROOT / "analysis/canonical/ppe/primary.pep.fa"],
    "gtr": [ROOT / "0903/karyotype/prepared/gtr/gtr.pep.fa"],
}


def read_fasta(path: Path):
    key = None
    seq = []
    for line in path.open(encoding="utf-8", errors="ignore"):
        if line.startswith(">"):
            if key:
                yield key, "".join(seq).replace("*", "")
            key = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if key:
        yield key, "".join(seq).replace("*", "")


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def parse_tree(text: str):
    text = text.strip().rstrip(";")
    pos = 0
    children = {}
    parent = {}
    ids = {}
    counter = [0]

    def node():
        i = counter[0]
        counter[0] += 1
        children[i] = []
        return i

    def parse():
        nonlocal pos
        n = node()
        if text[pos] == "(":
            pos += 1
            while True:
                c = parse()
                children[n].append(c)
                parent[c] = n
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            pos += 1
        name = text[start:pos].strip()
        if name:
            ids[n] = name
        return n

    root = parse()
    return root, children, parent, ids


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # ---- OrthoFinder OGs
    og_members = {}
    with OF.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = re.match(r"(OG\d+):\s*(.*)", line.strip())
            if not m:
                continue
            og_members[m.group(1)] = m.group(2).split()
    print(f"[turnover] OGs={len(og_members)}")

    # ---- candidate floral annotation -> frozen gIDs
    frozen_seq = {t: {} for t in TRACKS}
    tracker_for_species = {"pzi": ["pzi_A", "pzi_B"], "psp": ["psp_A", "psp_B"], "ppa": ["ppa_A", "ppa_B"],
                           "pcam": ["pcam"], "ppe": ["ppe"], "gtr": ["gtr"]}
    for sp, tracks in tracker_for_species.items():
        for t in tracks:
            p = FROZEN / f"{t}.fa"
            if not p.is_file():
                continue
            for gid, seq in read_fasta(p):
                frozen_seq[t][seq] = gid
    gene_og = {}
    for og, members in og_members.items():
        for m in members:
            gene_og[m] = og

    candidates = list(csv.DictReader(CAND.open(encoding="utf-8"), delimiter="\t"))
    og_families: dict[str, set[str]] = defaultdict(set)
    n_bridged = 0
    for c in candidates:
        if c.get("evidence_tier") not in ("tree_supported", "tree_weak"):
            continue
        sp = c.get("species")
        fam = c.get("subfamily_group") or c.get("family") or ""
        if not fam:
            continue
        # find the canonical peptide by record id
        rid = c.get("record_id", "")
        seq = None
        for path in CAND_PEP.get(sp, []):
            if not path.is_file():
                continue
            for gid, s in read_fasta(path):
                if gid == rid:
                    seq = s
                    break
            if seq:
                break
        if not seq:
            continue
        hit = None
        for t in tracker_for_species.get(sp, []):
            if seq in frozen_seq.get(t, {}):
                hit = frozen_seq[t][seq]
                break
        if not hit or hit not in gene_og:
            continue
        n_bridged += 1
        og_families[gene_og[hit]].add(fam)
    print(f"[turnover] candidate records bridged to OGs: {n_bridged}; floral OGs: {len(og_families)}")

    # ---- counts per species
    rows = []
    for og, members in og_members.items():
        counts = {sp: 0 for sp in SPECIES}
        per_track = {t: 0 for t in TRACKS}
        for m in members:
            t = m.split("__")[0]
            if t in per_track:
                per_track[t] += 1
        for sp, ts in SPECIES.items():
            counts[sp] = sum(per_track[t] for t in ts)
        n_species = sum(1 for sp in SPECIES if counts[sp] > 0)
        if n_species < 3:
            continue
        out_ref = [counts[s] for s in ("pcam", "ppe", "gtr", "rch")]
        out_med = float(np.median(out_ref)) if out_ref else 0.0
        trio = counts["pzi"] + counts["psp"] + counts["ppa"]
        trio_per = trio / 3.0
        expected = 2.0 * out_med
        ratio = trio_per / expected if expected > 0 else float("nan")
        rows.append({"og": og, **{f"n_{sp}": counts[sp] for sp in SPECIES},
                     "trio_total": trio, "outgroup_median": out_med, "retention_ratio": ratio,
                     "floral": ";".join(sorted(og_families.get(og, set())))})
    print(f"[turnover] OGs with >=3 species: {len(rows)}")

    # ---- tree-based least-squares ancestral states (equal branch lengths)
    root, children, parent, ids = parse_tree(TREE)
    order = sorted(children, key=lambda i: len(children[i]))
    # postorder
    post = []
    stack = [root]
    while stack:
        n = stack.pop()
        post.append(n)
        stack.extend(children[n])
    post = list(reversed(post))
    internal = [n for n in children if children[n]]
    leaf_ids = {n: ids[n] for n in children if not children[n] and n in ids}
    n_og = len(rows)
    vals = np.zeros((n_og, len(children)))
    for n, sp in leaf_ids.items():
        vals[:, n] = np.array([r[f"n_{sp}"] for r in rows], dtype=float)
    # iterative averaging on tree (equal branches)
    idx_map = {n: i for i, n in enumerate(children)}
    for _ in range(300):
        for n in internal:
            neigh = children[n] + ([parent[n]] if n in parent else [])
            vals[:, idx_map[n]] = vals[:, [idx_map[k] for k in neigh]].mean(axis=1)
    trio_node = None
    # find node whose leaf set == {pzi,psp,ppa}

    def leaves(n):
        if not children[n]:
            return {ids.get(n, "")}
        out = set()
        for c in children[n]:
            out |= leaves(c)
        return out

    for n in internal:
        if leaves(n) == {"pzi", "psp", "ppa"}:
            trio_node = n
            break
    stem_node = parent[trio_node] if trio_node is not None and trio_node in parent else None
    delta = None
    if trio_node is not None and stem_node is not None:
        delta = vals[:, idx_map[trio_node]] - vals[:, idx_map[stem_node]]
    for r, d in zip(rows, delta if delta is not None else np.full(n_og, np.nan)):
        r["delta_trio_stem"] = round(float(d), 4) if d == d else float("nan")
    write_tsv(OUT / "og_turnover.tsv", rows)

    # ---- floral enrichment
    floral_mask = np.array([bool(r["floral"]) for r in rows])
    ratio = np.array([r["retention_ratio"] for r in rows], dtype=float)
    delta_arr = np.array([r["delta_trio_stem"] for r in rows], dtype=float)
    ok = np.isfinite(ratio)
    print(f"[turnover] floral OGs tested: {int(floral_mask.sum())} of {len(rows)}")
    med_f = np.median(ratio[floral_mask & ok])
    med_b = np.median(ratio[(~floral_mask) & ok])
    u, p = mannwhitneyu(ratio[floral_mask & ok], ratio[(~floral_mask) & ok], alternative="two-sided")
    print(f"[turnover] retention ratio median floral={med_f:.3f} background={med_b:.3f} MWU p={p:.4g}")
    # expanded vs retained vs fractionated classes
    classes = np.where(ratio >= 1.3, "expanded", np.where(ratio < 0.7, "fractionated", "retained"))
    enr_rows = []
    for cls in ("expanded", "fractionated", "retained"):
        a = int(((classes == cls) & floral_mask).sum())
        b = int(((classes == cls) & ~floral_mask).sum())
        c = int((floral_mask).sum() - a)
        d = int((~floral_mask).sum() - b)
        orr, pf = fisher_exact([[a, b], [c, d]])
        enr_rows.append({"class": cls, "floral": a, "nonfloral": b,
                         "odds_ratio": round(orr, 4), "fisher_p": f"{pf:.4g}"})
    write_tsv(OUT / "turnover_floral_enrichment.tsv", enr_rows)
    print("[turnover] enrichment:", enr_rows)

    # top expanded floral OGs
    top = sorted([r for r in rows if r["floral"] and np.isfinite(r["retention_ratio"])],
                 key=lambda r: -r["retention_ratio"])[:20]
    write_tsv(OUT / "top_expanded_floral.tsv", top)
    for r in top[:10]:
        print("  ", r["og"], r["floral"], "ratio=", round(r["retention_ratio"], 2),
              "trio=", r["trio_total"], "out=", r["outgroup_median"], "delta=", r["delta_trio_stem"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
