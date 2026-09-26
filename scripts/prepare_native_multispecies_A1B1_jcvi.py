#!/usr/bin/env python3
"""Prepare native JCVI karyotype inputs for a PCAM-ordered multi-row figure."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SPECIES = ["pzi", "psp", "ppa", "pmal", "pav", "pcam", "pspe", "ppe", "pdu", "psal", "pmu_refseq", "par", "rch"]
POLY = {
    "pzi": ("pzi_hapA1", "pzi_hapB1"),
    "psp": ("psp_hapA1", "psp_hapB1"),
    "ppa": ("ppa_hapA", "ppa_hapB"),
}
GENERIC = ["pmal", "pav", "pcam", "pspe", "ppe", "pdu", "psal", "pmu_refseq", "par", "rch"]


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--order-table", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    return p.parse_args()


def read_order(path: Path) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["species"]].append(row)
    return out


def bed_records(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 6:
                    yield fields


def copy1_map(path: Path) -> dict[str, str]:
    out = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("copy1") and row.get("display_seqid"):
                out[row["copy1"]] = row["display_seqid"]
    return out


def rewrite_bed(records, seq_map: dict[str, str], out: Path) -> set[str]:
    out.parent.mkdir(parents=True, exist_ok=True)
    genes: set[str] = set()
    with out.open("w", encoding="utf-8") as handle:
        for fields in records:
            gene = fields[3]
            if gene in genes:
                continue
            fields[0] = seq_map.get(fields[0], fields[0])
            genes.add(gene)
            handle.write("\t".join(fields) + "\n")
    return genes


def read_simple(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) == 6:
                yield fields


def rewrite_simple(paths: list[Path], valid_q: set[str], valid_r: set[str], out: Path) -> int:
    seen: set[tuple[str, ...]] = set()
    count = 0
    with out.open("w", encoding="utf-8") as handle:
        for path in paths:
            if not path.is_file():
                continue
            for fields in read_simple(path):
                if fields[0] not in valid_q or fields[1] not in valid_q:
                    continue
                if fields[2] not in valid_r or fields[3] not in valid_r:
                    continue
                key = tuple(fields)
                if key in seen:
                    continue
                seen.add(key)
                handle.write("\t".join(fields) + "\n")
                count += 1
    return count


def main() -> int:
    a = args()
    inp = a.input_dir
    a.outdir.mkdir(parents=True, exist_ok=True)
    for sub in ("bed", "anchors"):
        (a.outdir / sub).mkdir(exist_ok=True)

    order = read_order(a.order_table)
    map_copy1 = copy1_map(inp / "maps" / "chromosome_rename_map.tsv")
    source_beds = {p.stem: p for p in (inp / "generic_bed").glob("*.bed")}
    source_beds.update({p.stem: p for p in (inp / "subgenome_bed").glob("*.bed")})

    valid: dict[str, set[str]] = {}
    seqids: dict[str, list[str]] = {}

    # PZI is a union BED: hapA1/B1 IDs support direct PZI-PSP/PPA edges,
    # while copy1 IDs support the existing PZI-copy1 external anchors.
    pzi_bed = a.outdir / "bed" / "pzi.bed"
    pzi_records = []
    for label in ("pzi_hapA1", "pzi_hapB1", "pzi_copy1"):
        for fields in bed_records(source_beds[label]):
            fields = list(fields)
            if label == "pzi_copy1":
                fields[0] = map_copy1.get(fields[0], fields[0])
            pzi_records.append(fields)
    valid["pzi"] = rewrite_bed(pzi_records, {}, pzi_bed)
    seqids["pzi"] = [r["display_label"] for r in sorted(order["pzi"], key=lambda x: int(x["display_order"]))]

    for species in ("psp", "ppa"):
        components = POLY[species]
        records = []
        for component in components:
            records.extend(bed_records(source_beds[component]))
        out = a.outdir / "bed" / f"{species}.bed"
        valid[species] = rewrite_bed(records, {}, out)
        seqids[species] = [r["display_label"] for r in sorted(order[species], key=lambda x: int(x["display_order"]))]

    for species in GENERIC:
        source = source_beds[species]
        mapping = {}
        # Generic chromosomes are shown as Chr01...Chr08 in their PCAM order.
        for row in order[species]:
            mapping[row["source_seqid"]] = row["display_label"]
        out = a.outdir / "bed" / f"{species}.bed"
        valid[species] = rewrite_bed(bed_records(source), mapping, out)
        seqids[species] = [r["display_label"] for r in sorted(order[species], key=lambda x: int(x["display_order"] ))]

    # Rewrite direct A/B anchors. Both files share the same PZI/target BED.
    anchor_sources = inp / "subgenome_anchors"
    direct = {
        "pzi_psp": [anchor_sources / "pzi_hapA1.psp_hapA1.anchors.simple", anchor_sources / "pzi_hapB1.psp_hapB1.anchors.simple"],
        "pzi_ppa": [anchor_sources / "pzi_hapA1.ppa_hapA.anchors.simple", anchor_sources / "pzi_hapB1.ppa_hapB.anchors.simple"],
    }
    anchor_counts = {}
    anchor_counts["pzi_psp"] = rewrite_simple(direct["pzi_psp"], valid["pzi"], valid["psp"], a.outdir / "anchors" / "pzi_psp.anchors.simple")
    anchor_counts["pzi_ppa"] = rewrite_simple(direct["pzi_ppa"], valid["pzi"], valid["ppa"], a.outdir / "anchors" / "pzi_ppa.anchors.simple")

    generic_anchor_dir = inp / "generic_anchors"
    for species in GENERIC:
        source = generic_anchor_dir / f"pzi_copy1.{species}.anchors.simple"
        anchor_counts[species] = rewrite_simple([source], valid["pzi"], valid[species], a.outdir / "anchors" / f"pzi_{species}.anchors.simple")

    # Native JCVI needs all sequence IDs in one line per layout row.
    (a.outdir / "multispecies_A1B1.seqids").write_text(
        "\n".join(",".join(seqids[s]) for s in SPECIES) + "\n", encoding="utf-8"
    )
    y_values = [0.96 - i * (0.88 / (len(SPECIES) - 1)) for i in range(len(SPECIES))]
    lines = ["# y, xstart, xend, rotation, color, label, va, bed"]
    for i, species in enumerate(SPECIES):
        label = species.replace("_refseq", "")
        lines.append(f" {y_values[i]:.5f}, .08, .92, 0, #4C5660, {label}, center, bed/{species}.bed")
    lines.append("# edges")
    edge_pairs = [("pzi", "psp", "pzi_psp"), ("pzi", "ppa", "pzi_ppa")]
    edge_pairs.extend(("pzi", species, species) for species in GENERIC)
    index = {s: i for i, s in enumerate(SPECIES)}
    for q, r, stem in edge_pairs:
        if anchor_counts.get(stem, 0):
            filename = stem if stem in {"pzi_psp", "pzi_ppa"} else f"pzi_{stem}"
            lines.append(f" e, {index[q]}, {index[r]}, anchors/{filename}.anchors.simple")
    (a.outdir / "multispecies_A1B1.layout").write_text("\n".join(lines) + "\n", encoding="utf-8")

    with (a.outdir / "native_input_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("edge\tanchor_file\tn_records\n")
        for stem, count in anchor_counts.items():
            filename = f"{stem}.anchors.simple" if stem in {"pzi_psp", "pzi_ppa"} else f"pzi_{stem}.anchors.simple"
            handle.write(f"pzi-{stem}\t{a.outdir / 'anchors' / filename}\t{count}\n")

    print(f"[OK] native JCVI input directory: {a.outdir}")
    print(f"[OK] rows: {len(SPECIES)}; PZI/PSP/PPA each contain 16 chromosomes")
    for stem, count in anchor_counts.items():
        print(f"[OK] {stem}: {count} anchor blocks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
