#!/usr/bin/env python3
"""Screen each target's top TE families for satellite-like misclassification.

Signals per family (consensus sequence):
* tandem periodicity: dominant inter-copy exact-k-mer offset in 100-2000 bp
  with >= 12 occurrences (satellites show strong periodic self-similarity)
* REXdb domain evidence: whether the family consensus hits REXdb proteins at
  all, and whether the hit family agrees with the EDTA class of the family

Outputs a TSV of flagged families and a per-class summary of how much bp sits
in suspected satellites.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("${PROJECT_ROOT}")
INPUTS = BASE / "0903/haploid_repartition/inputs"
OUT = BASE / "0903/figures/figure2_te_subfamily"
CACHE = Path("/tmp/te_subfamily")
TARGETS = {
    "PZI-A1": INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TEanno.gff3",
    "PZI-B1": INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TEanno.gff3",
    "PSP-A1": INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TEanno.gff3",
    "PSP-B1": INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TEanno.gff3",
    "PPA-A": INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TEanno.gff3",
    "PPA-B": INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TEanno.gff3",
    "PCAM": BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TEanno.gff3",
    "PPE": BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TEanno.gff3",
}
TELIBS = {
    "PZI-A1": INPUTS / "pzi/edta/hapA1/pzi_hapA1.fa.mod.EDTA.TElib.fa",
    "PZI-B1": INPUTS / "pzi/edta/hapB1/pzi_hapB1.fa.mod.EDTA.TElib.fa",
    "PSP-A1": INPUTS / "psp/edta/hapA1/psp_hapA1.fa.mod.EDTA.TElib.fa",
    "PSP-B1": INPUTS / "psp/edta/hapB1/psp_hapB1.fa.mod.EDTA.TElib.fa",
    "PPA-A": INPUTS / "ppa/edta/hapA/ppa_hapA.fa.mod.EDTA.TElib.fa",
    "PPA-B": INPUTS / "ppa/edta/hapB/ppa_hapB.fa.mod.EDTA.TElib.fa",
    "PCAM": BASE / "0903/TE/pcam/pcam.fa.mod.EDTA.TElib.fa",
    "PPE": BASE / "0903/TE/ppe/ppe.fa.mod.EDTA.TElib.fa",
}


def class_of(tok: str):
    t = tok.lower()
    if "gypsy" in t: return "Gypsy"
    if "copia" in t: return "Copia"
    if t.startswith("ltr/"): return "OtherLTR"
    if "helitron" in t: return "Helitron"
    if t.startswith("dna/") or t.startswith("mite/"): return "TIR/DNA"
    if t.startswith("line"): return "LINE"
    if t.startswith("sine"): return "SINE"
    if t == "unknown": return "OtherTE"
    return None


def read_fasta(path: Path) -> dict[str, str]:
    out, name, chunks = {}, None, []
    for line in path.open(encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks)
            name = line[1:].split("#")[0].split()[0]
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def periodicity(s: str, k: int = 15) -> tuple[int, int]:
    pos = {}
    for i in range(0, len(s) - k, 3):  # sample every 3rd position for speed
        pos.setdefault(s[i:i + k], []).append(i)
    hits = Counter()
    for ps in pos.values():
        if len(ps) > 60:
            continue
        for a in ps:
            for b in ps:
                if b > a:
                    hits[b - a] += 1
    best = 0
    cnt = 0
    for off, c in hits.most_common(200):
        if 100 <= off <= 2000 and c >= 12:
            best, cnt = off, c
            break
    return best, cnt


def main() -> int:
    rows = []
    for target, gff in TARGETS.items():
        fams = defaultdict(lambda: {"n": 0, "len": 0})
        with gff.open(encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                if raw.startswith("#") or not raw.strip():
                    continue
                f = raw.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                m = re.search(r"classification=([^;]+)", f[8])
                if not m:
                    continue
                cls = class_of(m.group(1).strip())
                if cls is None:
                    continue
                nm = re.search(r"Name=([^;]+)", f[8])
                fid = nm.group(1) if nm else "?"
                try:
                    s, e = int(f[3]), int(f[4])
                except ValueError:
                    continue
                d = fams[(cls, fid)]
                d["n"] += 1
                d["len"] += abs(e - s) + 1
        top = []
        by_cls = defaultdict(list)
        for (cls, fid), d in fams.items():
            by_cls[cls].append((d["len"], fid, d["n"]))
        for cls, items in by_cls.items():
            items.sort(reverse=True)
            for length, fid, n in items[:15]:
                if length >= 500_000:
                    top.append((cls, fid, length, n))

        lib = read_fasta(TELIBS[target])
        rex = {}
        cache = CACHE / f"{target}_rex.tsv"
        if cache.is_file():
            for r in csv.DictReader(cache.open(encoding="utf-8"), delimiter="\t"):
                rex[r["record_id"]] = (r["family"], r["sublineage"], r["clade"])

        for cls, fid, length, n in top:
            seq = lib.get(fid)
            if seq is None:
                base = re.sub(r"_(INT|LTR)$", "", fid)
                seq = lib.get(base)
            if not seq:
                continue
            per, cnt = periodicity(seq)
            rfam, rsub, rclade = rex.get(fid, ("", "", ""))
            flags = []
            if per and cnt >= 12:
                flags.append("tandem_period")
            if not rfam:
                flags.append("no_rexdb_domain")
            elif cls in ("LINE", "SINE", "OtherLTR", "Gypsy", "Copia") and rfam.startswith(("TIR", "Subclass")):
                flags.append("domain_class_mismatch")
            if flags:
                rows.append({
                    "target": target, "class": cls, "family": fid,
                    "raw_mb": round(length / 1e6, 2), "n_features": n,
                    "consensus_len": len(seq), "tandem_period": per, "period_count": cnt,
                    "rexdb_family": rfam, "rexdb_lineage": f"{rsub}/{rclade}",
                    "flags": ";".join(flags),
                })

    rows.sort(key=lambda r: (r["target"], r["class"], -r["raw_mb"]))
    fields = ["target", "class", "family", "raw_mb", "n_features", "consensus_len",
              "tandem_period", "period_count", "rexdb_family", "rexdb_lineage", "flags"]
    with (OUT / "suspected_misclassified_families.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # summary: suspected satellite bp per target/class
    agg = defaultdict(float)
    for r in rows:
        if "tandem_period" in r["flags"]:
            agg[(r["target"], r["class"])] += r["raw_mb"]
    print("== suspected satellite (tandem period) bp by target/class ==")
    for (t, c), mb in sorted(agg.items()):
        print(f"  {t:8s} {c:9s} {mb:7.1f} Mb (raw features; exclusive figure value is lower)")
    print(f"\nflagged families: {len(rows)} -> {OUT / 'suspected_misclassified_families.tsv'}")
    for r in rows[:20]:
        print("  ", r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
