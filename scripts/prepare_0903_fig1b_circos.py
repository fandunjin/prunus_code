#!/usr/bin/env python3
"""Prepare a five-species, Figure-1b-style Circos plot from 0903 only.

The script deliberately reads no legacy or non-0903 files.  The angular
coordinate system is the established PZI/PSP/PPA haplotype order.  PZI A-B
links occupy the central chord layer; the four adjacent JCVI species links
occupy the outer gaps between species rings.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ANGULAR_ORDER = (
    "C1A", "C6B", "C7B", "C2A", "C5B", "C3A", "C8B", "C7A",
    "C4A", "C1B", "C5A", "C2B", "C6A", "C3B", "C8A", "C4B",
)

# The diploid row order is read from the 0903 native JCVI seqids provenance.
# It maps accession-ordered PCAM/PPE chromosomes to the corresponding PZI A
# homologs while retaining the user's angular order for the Circos backbone.
DIPLOID_CANONICAL = ("C1A", "C3A", "C7A", "C8A", "C4A", "C5A", "C6A", "C2A")
SPECIES = ("pzi", "psp", "ppa", "pcam", "ppe")

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
RING_RADIUS = {
    "pzi": ("0.430r", "0.480r"),
    "psp": ("0.535r", "0.585r"),
    "ppa": ("0.640r", "0.690r"),
    "pcam": ("0.745r", "0.795r"),
    "ppe": ("0.850r", "0.900r"),
}
LINK_RADIUS = {
    "pzi_internal": ("0.350r", "0.060r"),
    "pzi_psp": ("0.505r", "0.455r"),
    "psp_ppa": ("0.610r", "0.560r"),
    "ppa_pcam": ("0.715r", "0.665r"),
    "pcam_ppe": ("0.820r", "0.770r"),
}


def parse_karyotype(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) < 6 or fields[0] != "chr":
            continue
        m = re.search(r"(C[1-8][AB])", fields[2])
        if m:
            lengths[m.group(1)] = int(fields[5])
    missing = [x for x in ANGULAR_ORDER if x not in lengths]
    if missing:
        raise ValueError(f"Missing PZI chromosomes in karyotype: {','.join(missing)}")
    return {x: lengths[x] for x in ANGULAR_ORDER}


def load_bed(path: Path) -> tuple[dict[str, tuple[str, int, int]], dict[str, int], list[str]]:
    genes: dict[str, tuple[str, int, int]] = {}
    max_end: dict[str, int] = defaultdict(int)
    seq_order: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[3] in genes:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            seqid = fields[0]
            if seqid not in max_end:
                seq_order.append(seqid)
            genes[fields[3]] = (seqid, start, end)
            max_end[seqid] = max(max_end[seqid], end)
    if not genes:
        raise ValueError(f"No genes found in BED: {path}")
    return genes, dict(max_end), seq_order


def parse_blocks(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("StartGeneA"):
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


def remap(value: int, source_length: int, target_length: int) -> int:
    return max(0, min(target_length, round(value * target_length / max(1, source_length))))


def canonical(seqid: str) -> str | None:
    match = re.search(r"(C[1-8][AB])", seqid)
    return match.group(1) if match else None


def build_diploid_map(seq_order: list[str]) -> dict[str, str]:
    if len(seq_order) < 8:
        raise ValueError(f"Expected eight diploid chromosomes, found {len(seq_order)}")
    return {seqid: DIPLOID_CANONICAL[i] for i, seqid in enumerate(seq_order[:8])}


def write_karyotype(data: Path, lengths: dict[str, int]) -> None:
    with (data / "karyotype.txt").open("w", encoding="ascii") as handle:
        for seqid in ANGULAR_ORDER:
            handle.write(
                f"chr - {seqid} {seqid} 0 {lengths[seqid]} {PAIR_COLORS[seqid[:2]]}\n"
            )


def write_species_bands(data: Path, lengths: dict[str, int], present: dict[str, set[str]]) -> None:
    for species in SPECIES:
        color = SPECIES_COLORS[species]
        with (data / f"{species}_band.txt").open("w", encoding="ascii") as handle:
            for seqid in ANGULAR_ORDER:
                if seqid in present[species]:
                    handle.write(
                        f"{seqid} 0 {lengths[seqid]} 1 fill_color={color},color={color}\n"
                    )


def write_internal_pzi(data: Path, lengths: dict[str, int], pzi_root: Path) -> int:
    source = pzi_root / "data" / "links_ab.txt"
    rows: list[str] = []
    pzi_lengths = parse_karyotype(pzi_root / "data" / "karyotype.txt")
    for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) < 6:
            continue
        q = canonical(fields[0])
        s = canonical(fields[3])
        if not q or not s or q not in lengths or s not in lengths:
            continue
        attrs = fields[6] if len(fields) > 6 else "color=pair1_a55,z=100,twist=no"
        rows.append(
            f"{q} {remap(int(fields[1]), pzi_lengths[q], lengths[q])} "
            f"{remap(int(fields[2]), pzi_lengths[q], lengths[q])} "
            f"{s} {remap(int(fields[4]), pzi_lengths[s], lengths[s])} "
            f"{remap(int(fields[5]), pzi_lengths[s], lengths[s])} {attrs}"
        )
    (data / "pzi_internal_links.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_adjacent_links(
    data: Path,
    native: Path,
    lengths: dict[str, int],
    species_pair: tuple[str, str],
    anchor_name: str,
    beds: dict[str, tuple[dict[str, tuple[str, int, int]], dict[str, int], list[str]]],
    seqmaps: dict[str, dict[str, str]],
    color_name: str,
) -> int:
    upper, lower = species_pair
    upper_genes, upper_ends, _ = beds[upper]
    lower_genes, lower_ends, _ = beds[lower]
    rows: list[str] = []
    for g1, g2, g3, g4, _score, orient in parse_blocks(native / "anchors" / anchor_name):
        if g1 not in upper_genes or g2 not in upper_genes or g3 not in lower_genes or g4 not in lower_genes:
            continue
        useq, us1, ue1 = upper_genes[g1]
        useq2, us2, ue2 = upper_genes[g2]
        lseq, ls1, le1 = lower_genes[g3]
        lseq2, ls2, le2 = lower_genes[g4]
        if useq != useq2 or lseq != lseq2:
            continue
        ucanon = seqmaps[upper].get(useq)
        lcanon = seqmaps[lower].get(lseq)
        if not ucanon or not lcanon or ucanon not in lengths or lcanon not in lengths:
            continue
        rows.append(
            f"{ucanon} {remap(min(us1, us2), upper_ends[useq], lengths[ucanon])} "
            f"{remap(max(ue1, ue2), upper_ends[useq], lengths[ucanon])} "
            f"{lcanon} {remap(min(ls1, ls2), lower_ends[lseq], lengths[lcanon])} "
            f"{remap(max(le1, le2), lower_ends[lseq], lengths[lcanon])} "
            f"color={color_name},z=50,twist={'yes' if orient == '-' else 'no'}"
        )
    destination = data / f"{upper}_{lower}_links.txt"
    destination.write_text("\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_config(out: Path) -> None:
    pair_colors = "\n".join(
        f"pair{i} = {color}" for i, color in enumerate(PAIR_COLORS.values(), 1)
    )
    species_colors = "\n".join(
        f"species_{name} = {color}" for name, color in SPECIES_COLORS.items()
    )
    plots = []
    for species in SPECIES:
        r0, r1 = RING_RADIUS[species]
        plots.append(
            f"<plot>\n type = histogram\n file = data/{species}_band.txt\n"
            f" r0 = {r0}\n r1 = {r1}\n min = 0\n max = 1\n"
            " orientation = out\n thickness = 0\n</plot>"
        )
    links = [
        "<link>\n file = data/pzi_internal_links.txt\n radius = 0.350r\n"
        " bezier_radius = 0.060r\n ribbon = yes\n thickness = 0.45p\n crest = 0.55\n</link>",
    ]
    link_specs = (
        ("pzi_psp", "0.505r", "0.455r"),
        ("psp_ppa", "0.610r", "0.560r"),
        ("ppa_pcam", "0.715r", "0.665r"),
        ("pcam_ppe", "0.820r", "0.770r"),
    )
    for name, radius, bezier in link_specs:
        links.append(
            f"<link>\n file = data/{name}_links.txt\n radius = {radius}\n"
            f" bezier_radius = {bezier}\n ribbon = yes\n thickness = 0.45p\n crest = 0.50\n</link>"
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
file = fig1b_pzi_psp_ppa_pcam_ppe_circos.png
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
    (out / "conf").mkdir(parents=True, exist_ok=True)
    (out / "conf" / "circos_fig1b.conf").write_text(config, encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("0903"))
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("0903/haploid_repartition/circos/fig1b_five_species"),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.out.resolve()
    native = root / "native_multispecies_A1B1_pcam_pairwise"
    pzi_circos = root / "haploid_repartition" / "circos" / "pzi"
    data = out / "data"
    data.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    lengths = parse_karyotype(pzi_circos / "data" / "karyotype.txt")

    beds = {
        species: load_bed(native / "bed" / f"{species}.bed")
        for species in SPECIES
    }
    seqmaps: dict[str, dict[str, str]] = {}
    present: dict[str, set[str]] = {}
    for species in ("pzi", "psp", "ppa"):
        seqmaps[species] = {seqid: seqid for seqid in beds[species][1] if seqid in lengths}
        present[species] = set(seqmaps[species].values())
    for species in ("pcam", "ppe"):
        seqmaps[species] = build_diploid_map(beds[species][2])
        present[species] = set(seqmaps[species].values())

    write_karyotype(data, lengths)
    write_species_bands(data, lengths, present)
    counts = {
        "pzi_internal": write_internal_pzi(data, lengths, pzi_circos),
        "pzi_psp": write_adjacent_links(data, native, lengths, ("pzi", "psp"), "pzi.psp.blocks.simple", beds, seqmaps, "species_psp_a45"),
        "psp_ppa": write_adjacent_links(data, native, lengths, ("psp", "ppa"), "psp.ppa.blocks.simple", beds, seqmaps, "species_ppa_a45"),
        "ppa_pcam": write_adjacent_links(data, native, lengths, ("ppa", "pcam"), "ppa.pcam.blocks.simple", beds, seqmaps, "species_pcam_a45"),
        "pcam_ppe": write_adjacent_links(data, native, lengths, ("pcam", "ppe"), "pcam.ppe.blocks.simple", beds, seqmaps, "species_ppe_a45"),
    }
    write_config(out)
    with (data / "fig1b_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("item", "value"))
        writer.writerow(("data_root", root))
        writer.writerow(("rings_inner_to_outer", ",".join(SPECIES)))
        writer.writerow(("angular_order", ",".join(ANGULAR_ORDER)))
        writer.writerow(("pzi_internal_source", pzi_circos / "data" / "links_ab.txt"))
        for name, count in counts.items():
            writer.writerow((f"{name}_links", count))
        for species in SPECIES:
            writer.writerow((f"{species}_chromosomes", ",".join(sorted(present[species]))))
    print(f"[OK] prepared {out}")
    for name, count in counts.items():
        print(f"[OK] {name}: {count} links")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
