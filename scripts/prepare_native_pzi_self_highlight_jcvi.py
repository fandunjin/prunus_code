#!/usr/bin/env python3
"""Prepare native JCVI karyotype inputs for a focused PZI self-synteny panel."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_native_multispecies_A1B1_pairwise_jcvi as native  # noqa: E402
import plot_pzi_self_highlight_multispecies as highlight  # noqa: E402


SPECIES = ["pzi_self_query", "pzi_self_target", "ppa", "psp", "ppe", "pcam", "gtr"]
EXTERNAL = ["ppa", "psp", "ppe", "pcam", "gtr"]
GRAY = "#A7B1BA"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=ROOT / "downstream_analysis/local_A1B1_ribbon_inputs")
    parser.add_argument("--self-bed", type=Path, default=ROOT / "downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/pzi_copy1.bed")
    parser.add_argument("--self-anchors", type=Path, default=ROOT / "downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/pzi_copy1.pzi_copy1.anchors.simple")
    parser.add_argument("--order-table", type=Path, default=ROOT / "analysis/karyotype/jcvi/figures/custom_svg/ribbons/multispecies_A1B1_pcam_aligned/chromosome_order_orientation.tsv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_pzi_self_highlight_ppa_psp_ppe_pcam_gtr")
    parser.add_argument("--min-score", type=int, default=0)
    parser.add_argument("--max-blocks-per-homolog", type=int, default=0)
    return parser.parse_args()


def map_pzi_order(order: dict[str, list[dict[str, str]]], copy_map: dict[str, str]) -> list[str]:
    rows = sorted(order["pzi"], key=lambda row: int(row["display_order"]))
    inverse = {display: copy1 for copy1, display in copy_map.items()}
    out = []
    for row in rows:
        copy1 = inverse[row["display_label"]]
        out.append(copy1 + ("-" if row["orientation_to_pcam"] == "reverse" else ""))
    return out


def stage_bed(input_path: Path, output_path: Path, seqid_map: dict[str, str] | None = None) -> None:
    rows = []
    with input_path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if seqid_map:
                fields[0] = seqid_map.get(fields[0], fields[0])
            rows.append("\t".join(fields))
    if not rows:
        raise ValueError(f"No BED records: {input_path}")
    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def read_order(path: Path) -> dict[str, list[dict[str, str]]]:
    rows: dict[str, list[dict[str, str]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.setdefault(row["species"], []).append(row)
    for species in rows:
        rows[species].sort(key=lambda row: int(row["display_order"]))
    return rows


def colored_self_blocks(self_anchors: Path, self_bed: Path, output: Path) -> list[dict[str, object]]:
    rows = highlight.read_highlights(self_anchors, self_bed)
    if not rows:
        raise ValueError("No selected PZI self blocks")
    with output.open("w", encoding="utf-8") as handle:
        for row in rows:
            qa, qb, ta, tb = row["genes"]
            handle.write(
                f'{row["color"]}*{qa}\t{qb}\t{ta}\t{tb}\t{row["score"]}\t{row["orientation"]}\n'
            )
    return rows


def gray_direct_blocks(
    input_dir: Path, self_bed: Path, target_bed: Path, filenames: list[str], output: Path,
) -> int:
    _, qgenes, _, _ = native.read_bed(self_bed)
    _, rgenes, _, _ = native.read_bed(target_bed)
    rows: list[str] = []
    for filename in filenames:
        path = input_dir / "subgenome_anchors" / filename
        for genes, score, orientation in native.parse_simple(path):
            if score < 0 or not all(gene in qgenes for gene in genes[:2]) or not all(gene in rgenes for gene in genes[2:]):
                continue
            rows.append(f'{GRAY}*' + "\t".join((*genes, str(score), orientation)))
    output.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return len(rows)


def projection_list(
    input_dir: Path, beds: dict, copy_map: dict[str, str], species: str,
) -> list[native.Proj]:
    specs = {
        "ppa": [
            ("pzi_hapA1", "ppa_hapA", "pzi_hapA1.ppa_hapA.anchors.simple"),
            ("pzi_hapB1", "ppa_hapB", "pzi_hapB1.ppa_hapB.anchors.simple"),
        ],
        "psp": [
            ("pzi_hapA1", "psp_hapA1", "pzi_hapA1.psp_hapA1.anchors.simple"),
            ("pzi_hapB1", "psp_hapB1", "pzi_hapB1.psp_hapB1.anchors.simple"),
        ],
        "ppe": [("pzi_copy1", "ppe", "pzi_copy1.ppe.anchors.simple")],
        "pcam": [("pzi_copy1", "pcam", "pzi_copy1.pcam.anchors.simple")],
        "gtr": [
            ("pzi_hapA1", "gtr", "pzi_hapA1.gtr.anchors.simple"),
            ("pzi_hapB1", "gtr", "pzi_hapB1.gtr.anchors.simple"),
        ],
    }
    out: list[native.Proj] = []
    for query, target, filename in specs[species]:
        folder = "generic_anchors" if query == "pzi_copy1" else "subgenome_anchors"
        out.extend(native.projection_blocks(
            query, species, target, input_dir / folder / filename,
            beds, copy_map, 0,
        ))
    return out


def gray_projection(upper: list[native.Proj], lower: list[native.Proj], output: Path, cap: int) -> int:
    temp = output.with_suffix(".tmp")
    count = native.write_pair_projection(upper, lower, temp, cap)
    lines = []
    if temp.exists():
        for raw in temp.read_text(encoding="utf-8").splitlines():
            lines.append(raw.split("*", 1)[1] if "*" in raw else raw)
    output.write_text("\n".join(f"{GRAY}*{line}" for line in lines) + ("\n" if lines else ""), encoding="utf-8")
    temp.unlink(missing_ok=True)
    return count


def write_audit(outdir: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "highlight_class", "source_anchor_line", "query_chr", "target_chr",
        "query_start_fraction", "query_end_fraction", "target_start_fraction",
        "target_end_fraction", "query_gene_a", "query_gene_b", "target_gene_a",
        "target_gene_b", "score", "orientation",
    ]
    with (outdir / "highlight_blocks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "highlight_class": row["class"], "source_anchor_line": row["line"],
                "query_chr": row["query_chr"], "target_chr": row["target_chr"],
                "query_start_fraction": f'{float(row["query_start"]):.6f}',
                "query_end_fraction": f'{float(row["query_end"]):.6f}',
                "target_start_fraction": f'{float(row["target_start"]):.6f}',
                "target_end_fraction": f'{float(row["target_end"]):.6f}',
                "query_gene_a": row["query_gene_a"], "query_gene_b": row["query_gene_b"],
                "target_gene_a": row["target_gene_a"], "target_gene_b": row["target_gene_b"],
                "score": row["score"], "orientation": row["orientation"],
            })


def main() -> int:
    args = parse_args()
    for path in (args.input_dir, args.self_bed, args.self_anchors, args.order_table):
        if not path.exists():
            raise FileNotFoundError(path)
    order = read_order(args.order_table)
    copy_map = native.copy1_map(args.input_dir / "maps" / "chromosome_rename_map.tsv")
    out = args.outdir
    (out / "bed").mkdir(parents=True, exist_ok=True)
    (out / "anchors").mkdir(parents=True, exist_ok=True)

    # Read existing JCVI BEDs and derive GTR ordering from the existing input.
    source_paths = {p.stem: p for p in (args.input_dir / "generic_bed").glob("*.bed")}
    source_paths.update({p.stem: p for p in (args.input_dir / "subgenome_bed").glob("*.bed")})
    beds = {label: native.read_bed(path) for label, path in source_paths.items()}
    native.derive_gtr_order(order, beds, args.input_dir)

    pzi_order = map_pzi_order(order, copy_map)
    stage_bed(args.self_bed, out / "bed/pzi_self_query.bed")
    stage_bed(args.self_bed, out / "bed/pzi_self_target.bed")

    components = {
        "ppa": ("ppa_hapA", "ppa_hapB"),
        "psp": ("psp_hapA1", "psp_hapB1"),
        "ppe": ("ppe",), "pcam": ("pcam",), "gtr": ("gtr",),
    }
    for species, labels in components.items():
        mapping = {row["source_seqid"]: row["display_label"] for row in order[species]}
        staged = []
        for label in labels:
            for fields in beds[label][0]:
                fields = list(fields)
                fields[0] = mapping.get(fields[0], fields[0])
                staged.append("\t".join(fields))
        (out / "bed" / f"{species}.bed").write_text("\n".join(staged) + "\n", encoding="utf-8")

    selected = colored_self_blocks(args.self_anchors, args.self_bed, out / "anchors/pzi_self_highlight.blocks.simple")
    gray_direct_blocks(
        args.input_dir, args.self_bed, out / "bed/ppa.bed",
        ["pzi_hapA1.ppa_hapA.anchors.simple", "pzi_hapB1.ppa_hapB.anchors.simple"],
        out / "anchors/pzi.ppa.blocks.simple",
    )
    projections = {species: projection_list(args.input_dir, beds, copy_map, species) for species in EXTERNAL}
    counts = [
        ("pzi-self", len(selected), "pzi_self_highlight.blocks.simple"),
        ("pzi-ppa", sum(1 for line in (out / "anchors/pzi.ppa.blocks.simple").read_text().splitlines() if line.strip()), "pzi.ppa.blocks.simple"),
    ]
    pair_order = [("ppa", "psp"), ("psp", "ppe"), ("ppe", "pcam"), ("pcam", "gtr")]
    for upper, lower in pair_order:
        path = out / "anchors" / f"{upper}.{lower}.blocks.simple"
        n = gray_projection(projections[upper], projections[lower], path, args.max_blocks_per_homolog)
        counts.append((f"{upper}-{lower}", n, path.name))

    seqids = [
        ",".join(pzi_order), ",".join(pzi_order),
    ]
    for species in EXTERNAL:
        seqids.append(",".join(
            row["display_label"] + ("-" if row["orientation_to_pcam"] == "reverse" else "")
            for row in order[species]
        ))
    (out / "pzi_self_highlight.seqids").write_text("\n".join(seqids) + "\n", encoding="utf-8")

    layout = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
    labels_beds = [
        ("PZI self query", "bed/pzi_self_query.bed"),
        ("PZI self target", "bed/pzi_self_target.bed"),
        ("ppa", "bed/ppa.bed"), ("psp", "bed/psp.bed"),
        ("ppe", "bed/ppe.bed"), ("pcam", "bed/pcam.bed"), ("gtr", "bed/gtr.bed"),
    ]
    ys = [0.96 - i * (0.86 / (len(labels_beds) - 1)) for i in range(len(labels_beds))]
    for y, (label, bed) in zip(ys, labels_beds):
        layout.append(f" {y:.5f}, .08, .92, 0, #4C5660, {label}, center, {bed}, center")
    layout.append("# edges")
    edge_files = [
        "pzi_self_highlight.blocks.simple", "pzi.ppa.blocks.simple",
        "ppa.psp.blocks.simple", "psp.ppe.blocks.simple", "ppe.pcam.blocks.simple",
        "pcam.gtr.blocks.simple",
    ]
    for i, filename in enumerate(edge_files):
        layout.append(f" e, {i}, {i+1}, anchors/{filename}")
    (out / "pzi_self_highlight.layout").write_text("\n".join(layout) + "\n", encoding="utf-8")
    with (out / "block_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("edge\tblock_count\tfile\n")
        for edge, count, filename in counts:
            handle.write(f"{edge}\t{count}\tanchors/{filename}\n")
    write_audit(out, selected)
    (out / "README.txt").write_text(
        "Native JCVI PZI self-highlight panel\n\n"
        "Rows: duplicated PZI copy1 reference tracks followed by ppa, psp, ppe, pcam and gtr.\n"
        "Colored edge: selected PZI self anchors; other edges are neutral gray JCVI blocks.\n"
        "PZI chromosome labels use copy1 chr01-chr16 and retain the existing PCAM order/orientation.\n",
        encoding="utf-8",
    )
    print(f"[OK] Prepared native JCVI inputs: {out}")
    print(f"[OK] highlighted PZI self blocks: {len(selected)}")
    for edge, count, filename in counts[1:]:
        print(f"[OK] {edge}: {count} blocks ({filename})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
