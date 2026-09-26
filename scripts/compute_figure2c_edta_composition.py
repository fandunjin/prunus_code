#!/usr/bin/env python3
"""Recompute per-haplotype TE composition from the 0903 per-haplotype EDTA GFFs.

Exclusive per-base assignment per chromosome: every covered base is assigned
once to the overlapping feature with the highest score density (score/bp),
then higher score, then longer span. Classes follow the Figure 2c legend.
"Other TEs" is the remainder of the total repeat union after the eight named
classes.
"""

import bisect
import csv
import re
from collections import defaultdict
from pathlib import Path

BASE = Path("${PROJECT_ROOT}")
INPUTS = BASE / "0903/haploid_repartition/inputs"
OUT = BASE / "0903/figures/figure2_abc"

TARGETS = {
    "pziA1": (INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TEanno.sum", "PZI-A1", "PZI-A"),
    "pziA2": (INPUTS / "pzi/edta/hapA2/pzi_hapA2.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapA2/pzi_hapA2.fa.mod.EDTA.TEanno.sum", "PZI-A2", "PZI-A"),
    "pziB1": (INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TEanno.sum", "PZI-B1", "PZI-B"),
    "pziB2": (INPUTS / "pzi/edta/hapB2/pzi_hapB2.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapB2/pzi_hapB2.fa.mod.EDTA.TEanno.sum", "PZI-B2", "PZI-B"),
    "pspA1": (INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TEanno.sum", "PSP-A1", "PSP-A"),
    "pspA2": (INPUTS / "psp/edta/hapA2/psp_hapA2.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapA2/psp_hapA2.fa.mod.EDTA.TEanno.sum", "PSP-A2", "PSP-A"),
    "pspB1": (INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TEanno.sum", "PSP-B1", "PSP-B"),
    "pspB2": (INPUTS / "psp/edta/hapB2/psp_hapB2.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapB2/psp_hapB2.fa.mod.EDTA.TEanno.sum", "PSP-B2", "PSP-B"),
    "ppaA": (INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TEanno.gff3", INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TEanno.sum", "PPA-A", "PPA-A"),
    "ppaB": (INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TEanno.gff3", INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TEanno.sum", "PPA-B", "PPA-B"),
    "pcam": (BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TEanno.gff3", BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TEanno.sum", "PCAM", "PCAM"),
    "ppe": (BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TEanno.gff3", BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TEanno.sum", "PPE", "PPE"),
}

NAMED = ["Copia", "Gypsy", "CACTA", "hAT", "Mutator", "PIF", "Tc1", "Helitron"]


def figure_class(cls):
    if cls == "LTR/Copia":
        return "Copia"
    if cls == "LTR/Gypsy":
        return "Gypsy"
    if cls in {"DNA/CMC-EnSpm", "DNA/DTC", "MITE/DTC"}:
        return "CACTA"
    if cls in {"DNA/DTA", "MITE/DTA", "DNA/hAT-Ac", "DNA/hAT-Tag1", "DNA/hAT-Tip100"}:
        return "hAT"
    if cls in {"DNA/DTM", "MITE/DTM", "DNA/MULE-MuDR"}:
        return "Mutator"
    if cls in {"DNA/DTH", "MITE/DTH", "DNA/PIF-Harbinger"}:
        return "PIF"
    if cls in {"DNA/DTT", "MITE/DTT"}:
        return "Tc1"
    if cls in {"DNA/Helitron", "RC/Helitron"}:
        return "Helitron"
    return "Other TEs"


def parse_class(attrs):
    m = re.search(r"classification=([^;]+)", attrs)
    return m.group(1).strip() if m else "unknown"


def genome_total(sum_path):
    with open(sum_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.match(r"Total Length:\s+(\d+)", line)
            if m:
                return int(m.group(1))
    raise ValueError(f"no Total Length in {sum_path}")


def union_bp(intervals):
    intervals = sorted(intervals)
    total = 0
    cs = ce = None
    for s, e in intervals:
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    if cs is not None:
        total += ce - cs + 1
    return total


def exclusive_by_chrom(feats_by_chrom):
    """Return (total_union, named_bp dict) using greedy priority assignment."""
    total_union = 0
    named = defaultdict(int)
    for chrom, feats in feats_by_chrom.items():
        total_union += union_bp([(s, e) for (s, e, _cls, _den, _sc, _ln) in feats])
        feats = sorted(feats, key=lambda f: (-f[3], -f[4], -f[5]))
        occupied = []
        starts = []

        def free_parts(s, e):
            parts = []
            i = bisect.bisect_right(starts, s) - 1
            if i < 0:
                i = 0
            cur = s
            while i < len(occupied) and occupied[i][0] <= e:
                os_, oe_ = occupied[i]
                if os_ > cur:
                    parts.append((cur, min(os_ - 1, e)))
                cur = max(cur, oe_ + 1)
                if cur > e:
                    break
                i += 1
            if cur <= e:
                parts.append((cur, e))
            return parts

        def insert(s, e):
            i = bisect.bisect_left(starts, s)
            occupied.insert(i, (s, e))
            starts.insert(i, s)
            j = i
            while j + 1 < len(occupied) and occupied[j + 1][0] <= occupied[j][1] + 1:
                occupied[j] = (occupied[j][0], max(occupied[j][1], occupied[j + 1][1]))
                del occupied[j + 1]
                del starts[j + 1]
            if j > 0 and occupied[j - 1][1] + 1 >= occupied[j][0]:
                occupied[j - 1] = (occupied[j - 1][0], max(occupied[j - 1][1], occupied[j][1]))
                del occupied[j]
                del starts[j]

        for s, e, cls, _den, _sc, _ln in feats:
            fc = figure_class(cls)
            if fc == "Other TEs":
                continue
            for ps, pe in free_parts(s, e):
                named[fc] += pe - ps + 1
                insert(ps, pe)
    return total_union, named


def main():
    per_hap = {}
    for target, (gff, sumf, disp, group) in TARGETS.items():
        genome = genome_total(sumf)
        feats_by_chrom = defaultdict(list)
        n = 0
        with open(gff, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9:
                    continue
                try:
                    s, e = int(f[3]), int(f[4])
                except ValueError:
                    continue
                if e < s:
                    s, e = e, s
                try:
                    score = float(f[5]) if f[5] not in (".", "") else 0.0
                except ValueError:
                    score = 0.0
                length = e - s + 1
                feats_by_chrom[f[0]].append((s, e, parse_class(f[8]), score / length, score, length))
                n += 1
        total_union, named = exclusive_by_chrom(feats_by_chrom)
        named_sum = sum(named.values())
        other = total_union - named_sum
        per_hap[target] = {"display": disp, "group": group, "genome": genome, "union": total_union,
                           "named": dict(named), "other": other}
        print(f"{target}: n={n} chroms={len(feats_by_chrom)} genome={genome/1e6:.1f}Mb "
              f"union={total_union/1e6:.1f}Mb ({100*total_union/genome:.2f}%) "
              f"named={named_sum/1e6:.1f}Mb other={other/1e6:.1f}Mb")

    fields = ["target", "display_label", "display_group", "class", "bp", "genome_bp", "genome_fraction"]
    with open(OUT / "figure2c_te_summary_edta.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for target, d in per_hap.items():
            for c in NAMED + ["Other TEs"]:
                bp = d["named"].get(c, 0) if c != "Other TEs" else d["other"]
                w.writerow({"target": target, "display_label": d["display"], "display_group": d["group"],
                            "class": c, "bp": bp, "genome_bp": d["genome"],
                            "genome_fraction": f"{bp/d['genome']:.8f}"})
    print("wrote", OUT / "figure2c_te_summary_edta.tsv")


if __name__ == "__main__":
    main()
