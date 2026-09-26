#!/usr/bin/env python3
"""Prepare a five-ring concentric Circos figure from the 0903 dataset.

Rings, inner to outer: PZI, PSP, PPA, PCAM, PPE.  Adjacent rings are connected
with the native JCVI blocks already stored under 0903, and the innermost PZI
ring also carries the PZI intra-subgenome A-B synteny links.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path


NATIVE = Path("0903/native_multispecies_A1B1_pcam_pairwise")
PZI_CIRCOS = Path("0903/haploid_repartition/circos/pzi")
DEFAULT_OUT = Path("0903/haploid_repartition/circos/concentric")

SPECIES = ("pzi", "psp", "ppa", "pcam", "ppe")

# Previous haploid Circos order: C1A ... C8B, grouped by homologous pair.
ANGULAR_ORDER = (
    "C1A", "C6B", "C7B", "C2A", "C5B", "C3A", "C8B", "C7A",
    "C4A", "C1B", "C5A", "C2B", "C6A", "C3B", "C8A", "C4B",
)

# PCAM/PPE chromosomes follow the diploid A-row order of the 0903 native
# pairwise layout; each is placed at its PZI canonical homologue.
DIPLOID_A_ORDER = ("C1A", "C3A", "C7A", "C8A", "C4A", "C5A", "C6A", "C2A")

RING_RADIUS = {
    "pzi": ("0.440r", "0.490r"),
    "psp": ("0.540r", "0.590r"),
    "ppa": ("0.640r", "0.690r"),
    "pcam": ("0.740r", "0.790r"),
    "ppe": ("0.840r", "0.890r"),
}

PAIR_COLORS = {
    "C1": "0,114,178",
    "C2": "230,159,0",
    "C3": "0,158,115",
    "C4": "213,94,0",
    "C5": "204,121,167",
    "C6": "86,180,233",
    "C7": "190,168,0",
    "C8": "80,80,80",
}

SPECIES_COLORS = {
    "pzi": "0,114,178",
    "psp": "0,158,115",
    "ppa": "230,159,0",
    "pcam": "213,94,0",
    "ppe": "204,121,167",
}

EDGES = (
    ("pzi", "psp", "pzi.psp.blocks.simple", "species_psp_a45", "0.520r", "0.502r"),
    ("psp", "ppa", "psp.ppa.blocks.simple", "species_ppa_a45", "0.620r", "0.602r"),
    ("ppa", "pcam", "ppa.pcam.blocks.simple", "species_pcam_a45", "0.720r", "0.702r"),
    ("pcam", "ppe", "pcam.ppe.blocks.simple", "species_ppe_a45", "0.820r", "0.802r"),
)


def parse_karyotype(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) >= 6 and fields[0] == "chr":
            match = re.search(r"(C[1-8][AB])", fields[2])
            if match and int(fields[5]) > lengths.get(match.group(1), -1):
                lengths[match.group(1)] = int(fields[5])
    return {chrom: lengths[chrom] for chrom in ANGULAR_ORDER if chrom in lengths}


def load_bed(path: Path):
    gene_locs: dict[str, tuple[str, int, int]] = {}
    max_end: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[3] in gene_locs:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            gene_locs[fields[3]] = (fields[0], start, end)
            max_end[fields[0]] = max(max_end[fields[0]], end)
    return gene_locs, max_end


def seqid_to_canonical(species: str, seqid: str) -> str | None:
    if species in ("pzi", "psp", "ppa"):
        return seqid if re.fullmatch(r"C[1-8][AB]", seqid) else None
    match = re.fullmatch(r"Chr(0?\d+)", seqid)
    if not match:
        return None
    index = int(match.group(1)) - 1
    if 0 <= index < len(DIPLOID_A_ORDER):
        return DIPLOID_A_ORDER[index]
    return None


def read_blocks(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6:
                continue
            try:
                yield fields[0], fields[1], fields[2], fields[3], int(float(fields[4])), fields[5]
            except ValueError:
                continue


def remap(value: int, source_len: int, target_len: int) -> int:
    return max(0, min(target_len, round(value * target_len / max(1, source_len))))


def build_species_links(
    out_data: Path,
    common: dict[str, int],
    native: Path,
) -> None:
    beds = {species: load_bed(native / "bed" / f"{species}.bed") for species in SPECIES}
    for upper, lower, filename, color, _radius, _bezier in EDGES:
        upper_locs, upper_ends = beds[upper]
        lower_locs, lower_ends = beds[lower]
        rows: list[str] = []
        skipped = 0
        for g1, g2, g3, g4, _score, orient in read_blocks(native / "anchors" / filename):
            if g1 not in upper_locs or g2 not in upper_locs or g3 not in lower_locs or g4 not in lower_locs:
                skipped += 1
                continue
            u1, u1s, u1e = upper_locs[g1]
            u2, u2s, u2e = upper_locs[g2]
            l1, l1s, l1e = lower_locs[g3]
            l2, l2s, l2e = lower_locs[g4]
            if u1 != u2 or l1 != l2:
                skipped += 1
                continue
            uchrom = seqid_to_canonical(upper, u1)
            lchrom = seqid_to_canonical(lower, l1)
            if uchrom not in common or lchrom not in common:
                skipped += 1
                continue
            ustart, uend = min(u1s, u2s), max(u1e, u2e)
            lstart, lend = min(l1s, l2s), max(l1e, l2e)
            ulen = max(1, upper_ends[u1])
            llen = max(1, lower_ends[l1])
            rows.append(
                f"{uchrom} {remap(ustart, ulen, common[uchrom])} "
                f"{remap(uend, ulen, common[uchrom])} "
                f"{lchrom} {remap(lstart, llen, common[lchrom])} "
                f"{remap(lend, llen, common[lchrom])} "
                f"color={color},z=50,twist={'yes' if orient == '-' else 'no'}"
            )
        (out_data / f"{upper}_{lower}_links.txt").write_text(
            "\n".join(rows) + "\n", encoding="ascii"
        )
        print(f"[OK] {upper}-{lower}: {len(rows)} links ({skipped} blocks skipped)")


def build_pzi_internal(
    out_data: Path,
    common: dict[str, int],
    pzi_root: Path,
) -> None:
    pzi_lengths = parse_karyotype(pzi_root / "data" / "karyotype.txt")
    rows: list[str] = []
    source = pzi_root / "data" / "links_ab.txt"
    if not source.is_file():
        raise FileNotFoundError(source)
    for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) < 7:
            continue
        qseq = re.search(r"(C[1-8][AB])", fields[0])
        sseq = re.search(r"(C[1-8][AB])", fields[3])
        if not qseq or not sseq:
            continue
        qc, sc = qseq.group(1), sseq.group(1)
        if qc not in common or sc not in common or qc not in pzi_lengths or sc not in pzi_lengths:
            continue
        qs, qe = int(fields[1]), int(fields[2])
        ss, se = int(fields[4]), int(fields[5])
        attrs = fields[6] if len(fields) > 6 else "color=pair1_a60,z=100,twist=no"
        rows.append(
            f"{qc} {remap(qs, pzi_lengths[qc], common[qc])} "
            f"{remap(qe, pzi_lengths[qc], common[qc])} "
            f"{sc} {remap(ss, pzi_lengths[sc], common[sc])} "
            f"{remap(se, pzi_lengths[sc], common[sc])} {attrs}"
        )
    (out_data / "pzi_internal_links.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
    print(f"[OK] pzi internal: {len(rows)} links")


def write_band_files(out_data: Path, common: dict[str, int]) -> None:
    for species in SPECIES:
        color = SPECIES_COLORS[species]
        with (out_data / f"{species}_bands.txt").open("w", encoding="ascii", newline="\n") as handle:
            for chrom in ANGULAR_ORDER:
                handle.write(
                    f"{chrom} 0 {common[chrom]} 1 "
                    f"fill_color={color},color={color}\n"
                )


def write_karyotype(out_data: Path, common: dict[str, int]) -> None:
    with (out_data / "karyotype.txt").open("w", encoding="ascii", newline="\n") as handle:
        for chrom in ANGULAR_ORDER:
            pair = chrom[:2]
            handle.write(
                f"chr - {chrom} {chrom} 0 {common[chrom]} {PAIR_COLORS[pair]}\n"
            )


def write_manifest(out_root: Path, common: dict[str, int]) -> None:
    lines = [
        "item\tvalue",
        "figure\tconcentric five-species CircOS",
        "rings_inner_to_outer\tpzi, psp, ppa, pcam, ppe",
        "angular_order\t" + ",".join(ANGULAR_ORDER),
        "intra_synteny\tPZI A-B links (pzi/data/links_ab.txt)",
        "inter_ring_synteny\tnative adjacent JCVI blocks under 0903/native_multispecies_A1B1_pcam_pairwise/anchors",
    ]
    for chrom in ANGULAR_ORDER:
        lines.append(f"{chrom}\t{common[chrom]}")
    (out_root / "data" / "concentric_manifest.tsv").write_text(
        "\n".join(lines) + "\n", encoding="ascii"
    )


def write_config(out_root: Path, common: dict[str, int]) -> None:
    species_colors = "\n".join(
        f"species_{species} = {color}" for species, color in SPECIES_COLORS.items()
    )
    pair_colors = "\n".join(
        f"pair{i} = {color}" for i, color in enumerate(PAIR_COLORS.values(), 1)
    )
    plots: list[str] = []
    for species in SPECIES:
        r0, r1 = RING_RADIUS[species]
        plots.append(
            f"<plot>\n type = histogram\n file = data/{species}_bands.txt\n "
            f"r0 = {r0}\n r1 = {r1}\n min = 0\n max = 1\n "
            "orientation = out\n thickness = 0\n</plot>"
        )
    links = [
        (
            "<link>\n file = data/pzi_internal_links.txt\n radius = 0.420r\n "
            "bezier_radius = 0.050r\n ribbon = yes\n thickness = 0.45p\n crest = 0.55\n</link>"
        )
    ]
    for upper, lower, _file, _color, radius, bezier in EDGES:
        links.append(
            f"<link>\n file = data/{upper}_{lower}_links.txt\n radius = {radius}\n "
            f"bezier_radius = {bezier}\n ribbon = yes\n thickness = 0.45p\n crest = 0.50\n</link>"
        )
    config = f"""<<include /etc/circos/colors_fonts_patterns.conf>>
<<include /etc/circos/housekeeping.conf>>
<colors>
{pair_colors}
{species_colors}
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {";".join(ANGULAR_ORDER)}
<image>
background = white
dir = figures
file = comparative_concentric_circos.png
png = yes
svg = yes
radius = 531.50p
angle_offset = -90
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.008r
</spacing>
radius = 0.930r
thickness = 1p
fill = no
stroke_thickness = 0
show_label = yes
label_radius = 1.045r
label_size = 13p
label_parallel = yes
</ideogram>
<plots>
{chr(10).join(plots)}
</plots>
<links>
{chr(10).join(links)}
</links>
"""
    (out_root / "conf").mkdir(parents=True, exist_ok=True)
    (out_root / "conf" / "circos_concentric.conf").write_text(config, encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-dir", type=Path, default=NATIVE)
    parser.add_argument("--pzi-circos-root", type=Path, default=PZI_CIRCOS)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out_root = args.out_root.resolve()
    native = args.native_dir.resolve()
    pzi_root = args.pzi_circos_root.resolve()
    data = out_root / "data"
    data.mkdir(parents=True, exist_ok=True)
    (out_root / "figures").mkdir(parents=True, exist_ok=True)

    common = parse_karyotype(pzi_root / "data" / "karyotype.txt")
    if len(common) != 16:
        raise ValueError(f"Expected 16 canonical chromosomes, got {len(common)}")
    write_karyotype(data, common)
    write_band_files(data, common)
    build_species_links(data, common, native)
    build_pzi_internal(data, common, pzi_root)
    write_config(out_root, common)
    write_manifest(out_root, common)
    print(f"[OK] concentric CircOS inputs: {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
