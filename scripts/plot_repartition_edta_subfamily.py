#!/usr/bin/env python3
"""Figure 2d-style TE composition subdivided to TE subfamily level.

NOTE (2026-09-17): the released figure2d / edta_composition_subfamily.tsv were
post-corrected to reclassify the mis-annotated LINE/L1 tandem satellite family
TE_00005699 as "Satellite — PziCent167"
(0903/fig1e_satdna/scripts/45_fig1e_fig2d_corrected.py). Re-running this script
reverts that correction; re-apply it afterwards.

Same per-base exclusive assignment as plot_repartition_edta_composition.py
(winner = highest RepeatMasker score, then span), but each interval is also
labelled with a subfamily:

* LTR/Gypsy and LTR/Copia families get a REXdb Viridiplantae v3.0 lineage
  (diamond blastx on the target TElib consensus sequences)
* TIR/DNA families are grouped by the EDTA class token
  (DTM/MULE -> Mutator; DTA/hAT -> hAT; DTC/CMC-EnSpm -> CACTA;
   DTH/PIF -> PIF/Harbinger; DTT/Tc1 -> Tc1/Mariner; else other TIR/DNA)
* LINE -> L1 / RTE / other; other classes stay single.

Outputs (0903/figures/figure2_te_subfamily/):
  edta_composition_subfamily.tsv   per target x class x subfamily bp
  figure2d_te_subfamily.{png,pdf,svg,tiff}
"""

from __future__ import annotations

import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

BASE = Path("${PROJECT_ROOT}")
REX = Path("/home/user/l1_validation/db/rexdb/repeatexplorer-rexdb-eeb8a7c")
DIAMOND = "${CONDA_ROOT}/bin/diamond"
OUT = BASE / "0903/figures/figure2_te_subfamily"
TMP = Path("/tmp/te_subfamily")
SID_RE = re.compile(r"^(?P<domain>.+?)__REXdb_ID(?P<num>\d+)$")

INPUTS = BASE / "0903/haploid_repartition/inputs"
TARGETS = [
    ("PZI-A1", INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TElib.fa"),
    ("PZI-B1", INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TEanno.gff3", INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TElib.fa"),
    ("PSP-A1", INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TElib.fa"),
    ("PSP-B1", INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TEanno.gff3", INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TElib.fa"),
    ("PPA-A", INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TEanno.gff3", INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TElib.fa"),
    ("PPA-B", INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TEanno.gff3", INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TElib.fa"),
    ("PCAM", BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TEanno.gff3", BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TElib.fa"),
    ("PPE", BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TEanno.gff3", BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TElib.fa"),
]

CLASSES = ["LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA", "Helitron", "LINE", "SINE", "Other TE"]
CLASS_COLORS = {
    "LTR/Gypsy": "#214E70", "LTR/Copia": "#743C5C", "Other LTR": "#9A6680",
    "TIR/DNA": "#176B6E", "Helitron": "#6F7F3B", "LINE": "#A65E2E",
    "SINE": "#C0834A", "Other TE": "#5D536F",
}
# subfamily shade ramps (dark -> light)
RAMPS = {
    "LTR/Gypsy": ["#214E70", "#3F7092", "#6392B0", "#8FB6CE", "#BBD3E3"],
    "LTR/Copia": ["#743C5C", "#95587A", "#B47A99", "#CDA2B8", "#E2CBD6"],
    "TIR/DNA": ["#176B6E", "#2F8C87", "#57ABA2", "#86C7BF", "#AEDACF", "#CFEAE6"],
    "LINE": ["#A65E2E", "#C17C4B", "#DBA882"],
    "Other LTR": ["#9A6680", "#B98BA3", "#D5BCC7"],
}
ORDER = {
    "LTR/Gypsy": ["OTA", "CRM", "Tekay", "other chromovirus", "other Gypsy"],
    "LTR/Copia": ["Ale", "SIRE", "Ivana", "Tork", "other Copia"],
    "TIR/DNA": ["Mutator", "hAT", "CACTA", "PIF/Harbinger", "Tc1/Mariner", "other TIR/DNA"],
    "LINE": ["L1", "RTE", "other LINE"],
    "Other LTR": ["Gypsy-like", "Copia-like", "unclassified"],
    "Helitron": ["Helitron"], "SINE": ["SINE"], "Other TE": ["Other TE"],
}


def genome_total(gff: Path) -> int:
    sum_path = gff.with_name(gff.name.replace(".gff3", ".sum"))
    with sum_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.match(r"Total Length:\s+(\d+)", line)
            if m:
                return int(m.group(1))
    raise ValueError(f"no Total Length in {sum_path}")


def normalize_class(token: str) -> str | None:
    """Verbatim port of normalize_category from plot_repartition_edta_composition.py
    (parent is always empty for GFF rows)."""
    t = token.lower().replace("_", "-")
    p = ""
    if "gypsy" in t:
        return "LTR/Gypsy"
    if "copia" in t:
        return "LTR/Copia"
    if t.startswith("ltr/") or p == "ltr":
        return "Other LTR"
    if "helitron" in t:
        return "Helitron"
    if t.startswith("tir/") or t.startswith("dna/") or p in {"tir", "dna"}:
        return "TIR/DNA"
    if t.startswith("line") or p == "line":
        return "LINE"
    if t.startswith("sine") or p == "sine" or "trna" in t:
        return "SINE"
    if t in {"snrna/unknown", "rrna/unknown", "repeat-fragment", "repeat_fragment"}:
        return "Other TE"
    if t.startswith("nonltr/") or t.startswith("non-tir/") or t.startswith("nonltr"):
        return "Other TE"
    if t.startswith("unknown/") or t == "unknown":
        return "Other TE"
    if "rdna" in t or t.startswith("45s") or t.startswith("rrna"):
        return None
    return "Other TE" if "/" in t else None


def dna_subfamily(token: str) -> str:
    t = token.lower()
    if "helitron" in t:
        return "Helitron"
    if any(k in t for k in ("dtm", "mule", "mutator")):
        return "Mutator"
    if any(k in t for k in ("dta", "hat")):
        return "hAT"
    if any(k in t for k in ("dtc", "cmc", "enspm", "cacta")):
        return "CACTA"
    if any(k in t for k in ("dth", "pif", "harbinger", "tourist")):
        return "PIF/Harbinger"
    if any(k in t for k in ("dtt", "tc1", "mariner")):
        return "Tc1/Mariner"
    return "other TIR/DNA"


def line_subfamily(token: str) -> str:
    t = token.lower()
    if t.startswith("line/l1"):
        return "L1"
    if "rte" in t:
        return "RTE"
    return "other LINE"


def rexdb_lineage(family_cls: str, info: tuple[str, ...]) -> str:
    """info = (class, order, family, sublineage, clade) from REXdb v3.0.

    For Ty3/gypsy rows the 6th column carry chromovirus/non-chromovirus and the
    clade (OTA, CRM, Tekay, ...). For Ty1/copia rows the lineage name
    (Ale/SIRE/Ivana/Tork/...) sits in the 4th column (sublineage)."""
    _cls, _order, fam, col4, col5 = info
    if family_cls == "LTR/Gypsy":
        clade = col5 or col4
        if clade == "OTA":
            return "OTA"
        if clade == "CRM":
            return "CRM"
        if clade == "Tekay":
            return "Tekay"
        if col4 == "chromovirus" and clade not in ("", "chromo-unclass", "chromo-outgroup"):
            return "other chromovirus"
        return "other Gypsy"
    if family_cls == "LTR/Copia":
        lineage = col4 or col5
        if lineage in ("Ale", "SIRE", "Ivana", "Tork"):
            return lineage
        return "other Copia"
    if family_cls == "Other LTR":
        if fam == "Ty3/gypsy":
            return "Gypsy-like"
        if fam == "Ty1/copia":
            return "Copia-like"
        return "unclassified"
    return ""


def load_rexdb_map() -> dict[str, tuple[str, ...]]:
    out = {}
    for line in (REX / "Viridiplantae_v3.0.classification").open(encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 5:
            vals = parts[1:6] + [""] * max(0, 5 - len(parts[1:6]))
            out[parts[0]] = tuple(vals)
    return out


def read_fasta(path: Path) -> dict[str, str]:
    out, name, chunks = {}, None, []
    for line in path.open(encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks)
            name = line[1:].split()[0].split("#")[0]
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def build_family_lineage(target: str, telib: Path, rex: dict) -> dict[str, str]:
    """family record name -> raw REXdb classification tuple keyed by record id."""
    cache = TMP / f"{target}_rex.tsv"
    if not cache.is_file():
        fa = TMP / f"{target}_telib.fa"
        seqs = read_fasta(telib)
        with fa.open("w", encoding="utf-8") as fh:
            for name, seq in seqs.items():
                fh.write(f">{name}\n{seq}\n")
        hits = TMP / f"{target}_hits.m8"
        subprocess.run([
            DIAMOND, "blastx", "-q", str(fa), "-d", str(REX / "Viridiplantae_v3.0.fasta.dmnd"),
            "-o", str(hits), "-f", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
            "--threads", "8", "--sensitive",
        ], check=True, capture_output=True)
        best: dict[str, tuple[float, str]] = {}
        for line in hits.open(encoding="utf-8", errors="replace"):
            q, s, pid, alen, ev, bs = line.rstrip("\n").split("\t")
            bs = float(bs)
            if q not in best or bs > best[q][0]:
                best[q] = (bs, s)
        with cache.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["record_id", "rexdb_id", "class", "order", "family", "sublineage", "clade"])
            for rid, (bs, s) in best.items():
                m = SID_RE.match(s)
                rexid = f"REXdb_ID{m.group('num')}" if m else ""
                info = rex.get(rexid, ("", "", "", "", ""))
                w.writerow([rid, rexid, *info])
    out = {}
    with cache.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["record_id"]] = (row["class"], row["order"], row["family"], row["sublineage"], row["clade"])
    return out


def parse_gff(gff: Path, lin: dict) -> dict[str, dict[str, dict[str, int]]]:
    """target -> seqid -> list handling; returns classes[class][subfamily] = bp."""
    classes: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    by_seq: dict[str, list] = defaultdict(list)
    attr_cls = re.compile(r"(?:^|;)classification=([^;]+)")
    attr_name = re.compile(r"(?:^|;)Name=([^;]+)")
    with gff.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                start, end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            if end <= start:
                continue
            m = attr_cls.search(f[8])
            token = m.group(1) if m else f[2]
            cls = normalize_class(token)
            if cls is None:
                continue
            fam_m = attr_name.search(f[8])
            fam = fam_m.group(1) if fam_m else ""
            try:
                score = float(f[5]) if f[5] not in {".", ""} else 0.0
            except ValueError:
                score = 0.0
            sub = subfamily_of(cls, token, fam, lin)
            by_seq[f[0]].append((start, end, score, end - start, cls, sub))

    for intervals in by_seq.values():
        points = sorted({p for s, e, *_ in intervals for p in (s, e)})
        starts, ends = defaultdict(list), defaultdict(list)
        for idx, (s, e, *_rest) in enumerate(intervals):
            starts[s].append(idx)
            ends[e].append(idx)
        active: set[int] = set()
        for left, right in zip(points, points[1:]):
            active.difference_update(ends[left])
            active.update(starts[left])
            if left == right or not active:
                continue
            winner = max(active, key=lambda i: (intervals[i][2], intervals[i][3]))
            cls, sub = intervals[winner][4], intervals[winner][5]
            classes[cls][sub] += right - left
    return classes


def subfamily_of(cls: str, token: str, fam: str, lin: dict) -> str:
    if cls in ("LTR/Gypsy", "LTR/Copia", "Other LTR"):
        info = lin.get(fam)
        if info is None and "_" in fam:
            base = re.sub(r"_(INT|LTR)$", "", fam)
            info = lin.get(base)
        if info:
            lab = rexdb_lineage(cls, info)
            if lab:
                return lab
        return "other Gypsy" if cls == "LTR/Gypsy" else ("other Copia" if cls == "LTR/Copia" else "unclassified")
    if cls == "TIR/DNA":
        return dna_subfamily(token)
    if cls == "LINE":
        return line_subfamily(token)
    return cls


def main() -> int:
    TMP.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    rex = load_rexdb_map()

    records = []
    for target, gff, telib in TARGETS:
        lin = build_family_lineage(target, telib, rex)
        classes = parse_gff(gff, lin)
        te = sum(v for sub in classes.values() for v in sub.values())
        genome = genome_total(gff)
        records.append({"target": target, "classes": classes, "te_bp": te, "genome_bp": genome})
        top = ", ".join(f"{c}={sum(classes[c].values())/1e6:.1f}Mb" for c in CLASSES if classes.get(c))
        print(f"{target}: TE={te/1e6:.1f}Mb  {top}")

    # TSV
    with (OUT / "edta_composition_subfamily.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["target", "class", "subfamily", "bp"])
        for rec in records:
            for cls in CLASSES:
                for sub, bp in sorted(rec["classes"].get(cls, {}).items(), key=lambda kv: -kv[1]):
                    w.writerow([rec["target"], cls, sub, bp])

    # plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7, "axes.linewidth": 0.75,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.facecolor": "white",
    })

    labels = [r["target"] for r in records]
    x = np.arange(len(labels))
    te_mb = np.array([r["te_bp"] / 1e6 for r in records])

    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    bottom = np.zeros(len(labels))
    legend_handles, legend_labels = [], []
    for cls in CLASSES:
        subs = ORDER.get(cls, [cls])
        present = [s for s in subs if any(r["classes"].get(cls, {}).get(s) for r in records)]
        extra = sorted({s for r in records for s in r["classes"].get(cls, {})} - set(present))
        ramp = RAMPS.get(cls, [CLASS_COLORS[cls]])
        for i, sub in enumerate(present + extra):
            vals = np.array([r["classes"].get(cls, {}).get(sub, 0) / 1e6 for r in records])
            color = ramp[i] if i < len(ramp) else ramp[-1]
            bars = ax.bar(x, vals, bottom=bottom, width=0.7, color=color, linewidth=0)
            bottom += vals
            legend_handles.append(bars[0])
            legend_labels.append(f"{cls} — {sub}" if sub != cls else cls)
    for i, rec in enumerate(records):
        pct = 100 * rec["te_bp"] / rec["genome_bp"]
        ax.text(i, bottom[i] + max(te_mb) * 0.015, f"{pct:.1f}%", ha="center", va="bottom", fontsize=5.6)
    ax.set_ylabel("TE-annotated sequence (Mb)")
    ax.set_ylim(0, max(te_mb) * 1.18)
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.tick_params(axis="both", length=2.5, width=0.6, pad=2, labelsize=6.5)
    ax.grid(axis="y", color="#E4E4E4", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(legend_handles, legend_labels, ncol=1, loc="center left",
              bbox_to_anchor=(1.01, 0.5), fontsize=5.0, columnspacing=0.8,
              handlelength=1.1, handleheight=0.8)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg", "tiff"):
        fig.savefig(OUT / f"figure2d_te_subfamily.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
