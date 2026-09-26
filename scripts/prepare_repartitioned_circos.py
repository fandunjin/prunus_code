#!/usr/bin/env python3
"""Prepare a revised PZI four-haplotype Circos figure from JCVI blocks."""

from __future__ import annotations

import argparse
from pathlib import Path

HAPS = ("A1", "A2", "B1", "B2")
ORDER = {
    "A1": [f"C{i}A" for i in range(1, 9)],
    "A2": [f"C{i}A" for i in range(8, 0, -1)],
    "B1": [f"C{i}B" for i in range(1, 9)],
    "B2": [f"C{i}B" for i in range(8, 0, -1)],
}
COLORS = {
    "C1A": "c1_ab", "C2A": "c2_ab", "C3A": "c3_ab", "C4A": "c4_ab",
    "C5A": "c5_ab", "C6A": "c6_ab", "C7A": "c7_ab", "C8A": "c8_ab",
}
RGB = {
    "c1_ab": "0,114,178", "c2_ab": "230,159,0", "c3_ab": "0,158,115",
    "c4_ab": "213,94,0", "c5_ab": "204,121,167", "c6_ab": "86,180,233",
    "c7_ab": "190,168,0", "c8_ab": "80,80,80",
}


def read_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    out = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) >= 4:
                out[fields[3]] = (fields[0], int(fields[1]), int(fields[2]))
    if not out:
        raise SystemExit(f"BED has no records: {path}")
    return out


def parse_blocks(path: Path, left_bed, right_bed):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 6:
                continue
            a, b, c, d, score, orient = fields
            if a not in left_bed or b not in left_bed or c not in right_bed or d not in right_bed:
                continue
            l1, l2 = left_bed[a], left_bed[b]
            r1, r2 = right_bed[c], right_bed[d]
            if l1[0] != l2[0] or r1[0] != r2[0]:
                continue
            yield (l1[0], min(l1[1], l2[1]), max(l1[2], l2[2]), r1[0], min(r1[1], r2[1]), max(r1[2], r2[2]), orient)


def genome_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    name = None
    length = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    sizes[name] = length
                name = raw[1:].split()[0]
                length = 0
            else:
                length += len(raw.strip())
    if name is not None:
        sizes[name] = length
    return sizes


def write_config(path: Path) -> None:
    color_lines = "\n".join(f"{name} = {rgb}" for name, rgb in RGB.items())
    chromosome_order = [f"{hap}_{chrom}" for hap in HAPS for chrom in ORDER[hap]]
    order_text = ";".join(chromosome_order)
    pair = f"{chromosome_order[0]} {chromosome_order[16]}"
    text = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
{color_lines}
within_gray = 165,165,165
within_gray_a48 = 165,165,165,a48
hapa = 55,95,145
hapb = 186,85,40
c1_ab_a65 = 0,114,178,a65
c2_ab_a65 = 230,159,0,a65
c3_ab_a65 = 0,158,115,a65
c4_ab_a65 = 213,94,0,a65
c5_ab_a65 = 204,121,167,a65
c6_ab_a65 = 86,180,233,a65
c7_ab_a65 = 190,168,0,a65
c8_ab_a65 = 80,80,80,a65
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {order_text}
<image>
background = white
dir = figures
file = pzi_repartition_circos.png
png = yes
svg = yes
radius = 1700p
angle_offset = -90
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.002r
<pairwise {pair}>
spacing = 10r
</pairwise>
</spacing>
radius = 0.900r
thickness = 26p
fill = yes
stroke_color = white
stroke_thickness = 1p
show_label = yes
label_radius = 1.035r
label_size = 15p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 1r
color = grey
thickness = 1p
multiplier = 1e-6
format = %d
<tick>
spacing = 10u
size = 6p
show_label = yes
label_size = 9p
label_offset = 2p
</tick>
</ticks>
<links>
<link>
file = data/links_background.txt
radius = 0.620r
bezier_radius = 0.070r
ribbon = yes
thickness = 1p
</link>
<link>
file = data/links_highlight.txt
radius = 0.610r
bezier_radius = 0.060r
ribbon = yes
thickness = 1p
</link>
</links>
"""
    path.write_text(text, encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--jcvi-root", type=Path, required=True)
    parser.add_argument("--pzi-input-root", type=Path, required=True)
    args = parser.parse_args()
    data = args.out_root / "data"
    conf = args.out_root / "conf"
    for directory in (data, conf, args.out_root / "figures"):
        directory.mkdir(parents=True, exist_ok=True)

    bed_cache = {hap: read_bed(args.input_root / f"pzi_{hap}.bed") for hap in HAPS}
    karyotype_rows = []
    for hap in HAPS:
        fasta = args.pzi_input_root / "pzi" / "genome" / f"pzi_hap{hap}.fa"
        for chrom, size in genome_sizes(fasta).items():
            karyotype_rows.append((f"{hap}_{chrom}", chrom, size, "hapa" if hap.startswith("A") else "hapb"))
    with (data / "karyotype.txt").open("w", encoding="ascii") as handle:
        for internal, label, size, color in karyotype_rows:
            handle.write(f"chr - {internal} {label} 0 {size} {color}\n")

    def write_links(filename: str, pairs, highlight: bool) -> int:
        count = 0
        with (data / filename).open("w", encoding="ascii") as handle:
            for left, right in pairs:
                simple = args.jcvi_root / f"pzi_{left}_pzi_{right}" / "blocks.simple"
                if not simple.is_file():
                    raise SystemExit(f"Missing JCVI block file: {simple}")
                for lchr, ls, le, rchr, rs, re, orient in parse_blocks(simple, bed_cache[left], bed_cache[right]):
                    color = "within_gray_a48"
                    if highlight:
                        color = f"{COLORS.get(lchr, 'c1_ab')}_a65"
                    twist = "yes" if orient == "-" else "no"
                    handle.write(f"{left}_{lchr} {ls} {le} {right}_{rchr} {rs} {re} color={color},z={'5' if not highlight else '100'},twist={twist}\n")
                    count += 1
        return count

    background_count = write_links("links_background.txt", (("A1", "A2"), ("B1", "B2")), False)
    highlight_count = write_links("links_highlight.txt", (("A1", "B1"), ("A2", "B2")), True)
    with (data / "link_summary.tsv").open("w", encoding="ascii") as handle:
        handle.write(f"layer\tblocks\nbackground\t{background_count}\nhighlight\t{highlight_count}\n")
    write_config(conf / "circos.conf")
    (data / "README.txt").write_text(
        "Internal Circos IDs include the phase (A1_/A2_/B1_/B2_) so repeated labels are unique.\n"
        "Ideogram labels remain C1A-C8A and C1B-C8B. Gray: within phase; color: A-B.\n",
        encoding="ascii",
    )
    print(f"[OK] Circos inputs: {args.out_root}; background={background_count}, highlight={highlight_count}")


if __name__ == "__main__":
    main()
