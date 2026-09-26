#!/usr/bin/env python3
"""Prepare one common-coordinate comparative CircOS overview for PZI/PSP/PPA."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SPECIES = ("pzi", "psp", "ppa")
HAPS = {
    "pzi": ("A1", "A2", "B1", "B2"),
    "psp": ("A1", "A2", "B1", "B2"),
    "ppa": ("A", "B"),
}
COMMON_ORDER = (
    "C1A", "C6B", "C7B", "C2A", "C5B", "C3A", "C8B", "C7A",
    "C4A", "C1B", "C5A", "C2B", "C6A", "C3B", "C8A", "C4B",
)
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
SPECIES_COLORS = {"pzi": "0,114,178", "psp": "0,158,115", "ppa": "213,94,0"}
BIN_SIZE = 1_000_000
EXTERNAL_ANCHOR_ROOT = Path(
    "downstream_analysis/local_A1B1_ribbon_inputs/subgenome_anchors"
)
EXTERNAL_BED_ROOT = Path("downstream_analysis/local_A1B1_ribbon_inputs/subgenome_bed")
WORK_ANCHOR_ROOT = Path(
    "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work"
)
WORK_BED_ROOT = Path("0903/downstream_analysis/figure2_repartition_ks/inputs")
PZI_CHROM_RENAME = Path(
    "downstream_analysis/local_A1B1_ribbon_inputs/maps/chromosome_rename_map.tsv"
)
PZI_REPARTITION_MAP = Path("0903/haploid_repartition/maps/pzi_repartition_map.tsv")
PSP_RENAME_MAP = Path("data/psp/haploid/psp_A1A2B1B2.rename_map.tsv")
PPA_RENAME_MAP = Path("data/ppa/haploid_rename/ppa_A_B.rename_map.tsv")


def parse_karyotype(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) >= 6 and fields[0] == "chr":
            result[fields[2]] = int(fields[5])
    if not result:
        raise ValueError(f"No chromosomes in {path}")
    return result


def canonical(seqid: str) -> str:
    match = re.search(r"(C[1-8][AB])", seqid)
    if not match:
        raise ValueError(f"Cannot parse canonical chromosome from {seqid}")
    return match.group(1)


def canonical_seqid(seqid: str) -> str:
    """Normalize local haplotype labels such as C01_A1 to C1A."""
    match = re.match(r"^C0?([1-8])_([AB])(?:[12])?$", seqid)
    if match:
        return f"C{match.group(1)}{match.group(2)}"
    return canonical(seqid) if re.search(r"C[1-8][AB]", seqid) else seqid


def karyotype_length_rows(circos_root: Path, species: str) -> list[tuple[str, int]]:
    """Return canonical chromosome labels and lengths for one species."""
    return [
        (canonical(seqid), length)
        for seqid, length in parse_karyotype(
            circos_root / species / "data" / "karyotype.txt"
        ).items()
    ]


def best_length_match(length: int, rows: list[tuple[str, int]]) -> tuple[str, int]:
    """Match an annotation span to the closest assembled chromosome length."""
    chrom, chrom_length = min(rows, key=lambda row: abs(row[1] - length))
    return chrom, chrom_length - length


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_pzi_old_to_final() -> dict[str, str]:
    """Map intermediate C0x_A1/B1 names to the final C1A-C8B system."""
    intermediate = read_tsv_rows(PZI_CHROM_RENAME)
    copy1_to_old: dict[str, str] = {}
    for row in intermediate:
        if row.get("copy1") and row.get("display_seqid"):
            copy1_to_old[row["copy1"]] = row["display_seqid"]
    final_rows = read_tsv_rows(PZI_REPARTITION_MAP)
    final_by_copy1 = {
        row["copy1_chr"]: row["new_seqid"] for row in final_rows if row.get("copy1_chr")
    }
    old_to_final = {
        old: final_by_copy1[copy1]
        for copy1, old in copy1_to_old.items()
        if copy1 in final_by_copy1
    }
    return old_to_final


def build_local_to_final() -> dict[str, dict[str, str]]:
    """Map the local A1B1 bundle seqids to the final canonical labels."""
    pzi_old_to_final = build_pzi_old_to_final()
    psp_map: dict[str, str] = {}
    for row in read_tsv_rows(PSP_RENAME_MAP):
        if row.get("new_seqid") and row.get("pzi_target") in pzi_old_to_final:
            psp_map[row["new_seqid"]] = pzi_old_to_final[row["pzi_target"]]
    ppa_map: dict[str, str] = {}
    for row in read_tsv_rows(PPA_RENAME_MAP):
        label = row.get("pzi_copy1_haplotype_label")
        if row.get("renamed_seqid") and label in pzi_old_to_final:
            ppa_map[row["renamed_seqid"]] = pzi_old_to_final[label]
    pzi_map = {
        old: final for old, final in pzi_old_to_final.items()
    }
    return {"pzi": pzi_map, "psp": psp_map, "ppa": ppa_map}


def pair_color(seqid: str) -> str:
    return PAIR_COLORS[canonical(seqid)[:2]]


def all_lengths(circos_root: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for species in SPECIES:
        for seqid, length in parse_karyotype(circos_root / species / "data" / "karyotype.txt").items():
            key = canonical(seqid)
            lengths[key] = max(lengths.get(key, 0), length)
    return lengths


def source_seqids(circos_root: Path, species: str) -> dict[str, tuple[str, int]]:
    result: dict[str, tuple[str, int]] = {}
    for seqid, length in parse_karyotype(circos_root / species / "data" / "karyotype.txt").items():
        result[seqid] = (canonical(seqid), length)
    return result


def remap(value: int, source_len: int, target_len: int) -> int:
    return max(0, min(target_len, round(value * target_len / source_len)))


def remap_tracks(circos_root: Path, out_data: Path, common_lengths: dict[str, int]) -> None:
    out_data.mkdir(parents=True, exist_ok=True)
    track_names = ("gene_density", "te_coverage", "ncrna_density", "gc_fraction")
    for species in SPECIES:
        seqmap = source_seqids(circos_root, species)
        for track in track_names:
            values: dict[str, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
            source = circos_root / species / "data" / f"{track}.txt"
            if not source.is_file():
                raise FileNotFoundError(source)
            for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
                fields = raw.split()
                if len(fields) < 4 or fields[0] not in seqmap:
                    continue
                seqid, start, end, value = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
                chrom, source_len = seqmap[seqid]
                target_len = common_lengths[chrom]
                midpoint = remap((start + end) // 2, source_len, target_len)
                target_bin = (midpoint // BIN_SIZE) * BIN_SIZE
                values[chrom][target_bin].append(value)
            destination = out_data / f"{species}_{track}.txt"
            with destination.open("w", encoding="ascii", newline="\n") as handle:
                for chrom in COMMON_ORDER:
                    length = common_lengths[chrom]
                    for start in range(0, length, BIN_SIZE):
                        end = min(start + BIN_SIZE, length)
                        row = values[chrom].get(start, [])
                        value = sum(row) / len(row) if row else 0.0
                        handle.write(f"{chrom} {start} {end} {value:.6f}\n")


def remap_links(circos_root: Path, out_data: Path, common_lengths: dict[str, int]) -> None:
    for species in SPECIES:
        out_path = out_data / f"{species}_links.txt"
        rows: list[str] = []
        seqmap = source_seqids(circos_root, species)
        source = circos_root / species / "data" / "links_ab.txt"
        for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = raw.split()
            if len(fields) < 7 or fields[0] not in seqmap or fields[3] not in seqmap:
                continue
            qid, qstart, qend, sid, sstart, send = fields[:6]
            qchrom, qlen = seqmap[qid]
            schrom, slen = seqmap[sid]
            attrs = fields[6]
            attrs = re.sub(r"color=[^,]+", f"color=species_{species}_a55", attrs)
            rows.append(
                f"{qchrom} {remap(int(qstart), qlen, common_lengths[qchrom])} "
                f"{remap(int(qend), qlen, common_lengths[qchrom])} "
                f"{schrom} {remap(int(sstart), slen, common_lengths[schrom])} "
                f"{remap(int(send), slen, common_lengths[schrom])} {attrs}"
            )
        out_path.write_text("\n".join(rows) + "\n", encoding="ascii")


def read_gene_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    """Read JCVI BED coordinates keyed by gene ID."""
    genes: dict[str, tuple[str, int, int]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[3] in genes:
                continue
            try:
                genes[fields[3]] = (fields[0], int(fields[1]), int(fields[2]))
            except ValueError:
                continue
    if not genes:
        raise ValueError(f"No genes found in BED: {path}")
    return genes


def strip_gene_prefix(gene: str, prefixes: tuple[str, ...]) -> str:
    for prefix in prefixes:
        if gene.startswith(prefix):
            return gene[len(prefix):]
    return gene


def normalize_gene_bed(
    path: Path, prefixes: tuple[str, ...]
) -> dict[str, tuple[str, int, int]]:
    return {
        strip_gene_prefix(gene, prefixes): record
        for gene, record in read_gene_bed(path).items()
    }


def parse_anchor_rows(path: Path):
    """Yield query/subject gene pairs, score and orientation from JCVI simple anchors."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6:
                continue
            try:
                yield fields[0], fields[1], fields[2], fields[3], int(float(fields[4])), fields[5]
            except ValueError:
                continue


def remap_external_links(
    out_data: Path,
    common_lengths: dict[str, int],
    circos_root: Path,
    local_anchor_root: Path,
    local_bed_root: Path,
    work_anchor_root: Path,
    work_bed_root: Path,
) -> None:
    """Create PZI-to-reference links on the common C1A-C8B angular system.

    PCAM/PPE chromosomes are not added as a second karyotype.  Their anchor
    coordinates are therefore projected to the dominant PZI canonical
    chromosome supported by the same reference chromosome.  PSP/PPA use the
    direct local JCVI bundle and are matched to the final canonical naming by
    chromosome length, because that bundle predates the final renaming.
    """
    specs = {
        "psp": {
            "source": "local",
            "local_to_current": True,
            "queries": (
                ("pzi_hapA1", "pzi_hapA1.psp_hapA1.anchors.simple"),
                ("pzi_hapB1", "pzi_hapB1.psp_hapB1.anchors.simple"),
            ),
            "reference_beds": ("psp_hapA1", "psp_hapB1"),
            "reference_prefixes": ("psp_hapA1_", "psp_hapB1_"),
            "color": "species_psp_a55",
        },
        "ppa": {
            "source": "local",
            "local_to_current": True,
            "queries": (
                ("pzi_hapA1", "pzi_hapA1.ppa_hapA.anchors.simple"),
                ("pzi_hapB1", "pzi_hapB1.ppa_hapB.anchors.simple"),
            ),
            "reference_beds": ("ppa_hapA", "ppa_hapB"),
            "reference_prefixes": (),
            "color": "species_ppa_a55",
        },
        "pcam": {
            "source": "work",
            "local_to_current": False,
            "queries": (
                ("pzi_hapA1", "pzi_hapA1.pcam.anchors.simple"),
                ("pzi_hapB1", "pzi_hapB1.pcam.anchors.simple"),
            ),
            "reference_beds": ("pcam",),
            "reference_prefixes": (),
            "color": "species_pcam_a55",
        },
        "ppe": {
            "source": "work",
            "local_to_current": False,
            "queries": (
                ("pzi_hapA1", "pzi_hapA1.ppe.anchors.simple"),
                ("pzi_hapB1", "pzi_hapB1.ppe.anchors.simple"),
            ),
            "reference_beds": ("ppe",),
            "reference_prefixes": (),
            "color": "species_ppe_a55",
        },
    }
    query_prefixes = ("pzi_hapA1_", "pzi_hapB1_")
    local_to_final = build_local_to_final()
    audits: list[tuple[str, str, str, int, int, str, str]] = []
    for species, spec in specs.items():
        query_specs = spec["queries"]
        color = spec["color"]
        if spec["source"] == "local":
            anchor_root = local_anchor_root
            bed_root = local_bed_root
            reference_paths = [bed_root / f"{label}.bed" for label in spec["reference_beds"]]
        else:
            anchor_root = work_anchor_root
            bed_root = work_bed_root
            reference_paths = [
                bed_root / label / f"{label}.bed" for label in spec["reference_beds"]
            ]
        target_bed: dict[str, tuple[str, int, int]] = {}
        for path in reference_paths:
            for gene, record in normalize_gene_bed(path, spec["reference_prefixes"]).items():
                target_bed.setdefault(gene, record)
        query_beds = {}
        for label, _ in query_specs:
            if spec["source"] == "local":
                bed_path = bed_root / f"{label}.bed"
            else:
                bed_path = bed_root / label / f"{label}.bed"
            query_beds[label] = normalize_gene_bed(bed_path, query_prefixes)
        source_rows: list[tuple[str, str, str, str, int, str, dict[str, tuple[str, int, int]], dict[str, tuple[str, int, int]]]] = []
        target_votes: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for label, filename in query_specs:
            anchor_path = anchor_root / f"{label}_vs_{species}" / filename
            if not anchor_path.is_file():
                anchor_path = anchor_root / filename
            if not anchor_path.is_file():
                raise FileNotFoundError(anchor_path)
            qbed = query_beds[label]
            for raw_q1, raw_q2, raw_t1, raw_t2, score, orient in parse_anchor_rows(anchor_path):
                q1 = strip_gene_prefix(raw_q1, query_prefixes)
                q2 = strip_gene_prefix(raw_q2, query_prefixes)
                t1 = strip_gene_prefix(raw_t1, spec["reference_prefixes"])
                t2 = strip_gene_prefix(raw_t2, spec["reference_prefixes"])
                if q1 not in qbed or q2 not in qbed or t1 not in target_bed or t2 not in target_bed:
                    continue
                qseq, _, _ = qbed[q1]
                qseq2, _, _ = qbed[q2]
                tseq, _, _ = target_bed[t1]
                tseq2, _, _ = target_bed[t2]
                qseq = canonical_seqid(qseq)
                qseq2 = canonical_seqid(qseq2)
                tseq = canonical_seqid(tseq)
                tseq2 = canonical_seqid(tseq2)
                if qseq != qseq2 or tseq != tseq2 or qseq not in common_lengths:
                    continue
                target_votes[tseq][qseq] += max(1, score)
                source_rows.append((q1, q2, t1, t2, score, orient, qbed, target_bed))
        target_map = {
            tseq: max(votes.items(), key=lambda item: item[1])[0]
            for tseq, votes in target_votes.items()
            if votes
        }
        query_lengths = {
            canonical_seqid(seqid): max(
                value[2]
                for bed in query_beds.values()
                for value in bed.values()
                if canonical_seqid(value[0]) == canonical_seqid(seqid)
            )
            for seqid in {value[0] for bed in query_beds.values() for value in bed.values()}
        }
        target_lengths = {
            canonical_seqid(seqid): max(
                value[2]
                for value in target_bed.values()
                if canonical_seqid(value[0]) == canonical_seqid(seqid)
            )
            for seqid in {value[0] for value in target_bed.values()}
        }
        rows: list[str] = []
        for q1, q2, t1, t2, score, orient, qbed, tbed in source_rows:
            raw_qseq, qs1, qe1 = qbed[q1]
            _qseq2, qs2, qe2 = qbed[q2]
            raw_tseq, ts1, te1 = tbed[t1]
            _tseq2, ts2, te2 = tbed[t2]
            qseq = canonical_seqid(raw_qseq)
            _qseq2 = canonical_seqid(_qseq2)
            tseq = canonical_seqid(raw_tseq)
            _tseq2 = canonical_seqid(_tseq2)
            if spec["local_to_current"]:
                qseq_axis = local_to_final["pzi"].get(raw_qseq)
                target_chrom = local_to_final[species].get(raw_tseq)
            else:
                qseq_axis = qseq
                target_chrom = target_map.get(tseq)
            if not target_chrom:
                continue
            qstart, qend = min(qs1, qs2), max(qe1, qe2)
            tstart, tend = min(ts1, ts2), max(te1, te2)
            target_max = max(1, target_lengths[tseq])
            rows.append(
                f"{qseq_axis} {remap(qstart, max(1, query_lengths[qseq]), common_lengths[qseq_axis])} "
                f"{remap(qend, max(1, query_lengths[qseq]), common_lengths[qseq_axis])} "
                f"{target_chrom} {remap(tstart, target_max, common_lengths[target_chrom])} "
                f"{remap(tend, target_max, common_lengths[target_chrom])} "
                f"color={color},z=20,twist={'yes' if orient == '-' else 'no'}"
            )
        (out_data / f"pzi_{species}_links.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
        if spec["local_to_current"]:
            for local_seqid, target_chrom in sorted(local_to_final[species].items()):
                audits.append(
                    (
                        species,
                        local_seqid,
                        target_chrom,
                        0,
                        0,
                        ";".join(path.name for path in reference_paths),
                        "local_rename_map",
                    )
                )
        else:
            for tseq, votes in sorted(target_votes.items()):
                target_chrom = target_map.get(tseq, "")
                audits.append(
                    (
                        species,
                        tseq,
                        target_chrom,
                        len(votes),
                        sum(votes.values()),
                        ";".join(path.name for path in reference_paths),
                        "vote_based_projection",
                    )
                )
    with (out_data / "external_reference_mapping.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "reference",
                "reference_seqid",
                "projected_canonical",
                "candidate_count",
                "weighted_anchor_score_or_length_delta",
                "bed",
                "projection",
            )
        )
        writer.writerows(audits)


def write_karyotype(out_data: Path, common_lengths: dict[str, int]) -> None:
    with (out_data / "karyotype.txt").open("w", encoding="ascii", newline="\n") as handle:
        for chrom in COMMON_ORDER:
            pair = chrom[:2]
            handle.write(f"chr - {chrom} {chrom} 0 {common_lengths[chrom]} {PAIR_COLORS[pair]}\n")


def write_config(out_root: Path, common_lengths: dict[str, int]) -> None:
    plot_specs = (
        ("pzi", "gene_density", "heatmap", "0.945r", "0.985r", "0", "205", "color = spectral-11-div"),
        ("pzi", "te_coverage", "histogram", "0.900r", "0.938r", "0", "1", "fill_color = 112,135,52\\ncolor = 112,135,52"),
        ("pzi", "ncrna_density", "histogram", "0.855r", "0.893r", "0", "1000", "fill_color = 0,114,178\\ncolor = 0,114,178"),
        ("pzi", "gc_fraction", "line", "0.835r", "0.850r", "0.34", "0.56", "color = 42,88,108"),
        ("psp", "gene_density", "heatmap", "0.780r", "0.820r", "0", "205", "color = spectral-11-div"),
        ("psp", "te_coverage", "histogram", "0.735r", "0.773r", "0", "1", "fill_color = 112,135,52\\ncolor = 112,135,52"),
        ("psp", "ncrna_density", "histogram", "0.690r", "0.728r", "0", "1000", "fill_color = 0,158,115\\ncolor = 0,158,115"),
        ("psp", "gc_fraction", "line", "0.670r", "0.685r", "0.34", "0.56", "color = 42,88,108"),
        ("ppa", "gene_density", "heatmap", "0.615r", "0.655r", "0", "205", "color = spectral-11-div"),
        ("ppa", "te_coverage", "histogram", "0.570r", "0.608r", "0", "1", "fill_color = 112,135,52\\ncolor = 112,135,52"),
        ("ppa", "ncrna_density", "histogram", "0.525r", "0.563r", "0", "1000", "fill_color = 213,94,0\\ncolor = 213,94,0"),
        ("ppa", "gc_fraction", "line", "0.505r", "0.520r", "0.34", "0.56", "color = 42,88,108"),
    )
    plots: list[str] = []
    for species, track, kind, r0, r1, minimum, maximum, options in plot_specs:
        plots.append(
            f"<plot>\n type = {kind}\n file = data/{species}_{track}.txt\n r0 = {r0}\n r1 = {r1}\n min = {minimum}\n max = {maximum}\n {options.replace('\\n', chr(10) + ' ')}\n thickness = 0\n</plot>"
            if kind != "line"
            else
            f"<plot>\n type = line\n file = data/{species}_{track}.txt\n r0 = {r0}\n r1 = {r1}\n min = {minimum}\n max = {maximum}\n {options}\n thickness = 1p\n</plot>"
        )
    link_specs = (
        ("pzi_links.txt", "0.490r", "0.410r", "0.55p", "0.42"),
        ("pzi_psp_links.txt", "0.402r", "0.322r", "0.50p", "0.40"),
        ("pzi_ppa_links.txt", "0.314r", "0.234r", "0.50p", "0.40"),
        ("pzi_pcam_links.txt", "0.226r", "0.146r", "0.48p", "0.38"),
        ("pzi_ppe_links.txt", "0.138r", "0.058r", "0.48p", "0.38"),
    )
    links = "\n".join(
        (
            f"<link>\n file = data/{filename}\n radius = {radius}\n "
            f"bezier_radius = {bezier}\n ribbon = yes\n thickness = {thickness}\n "
            f"crest = {crest}\n</link>"
        )
        for filename, radius, bezier, thickness, crest in link_specs
    )
    colors = "\n".join(f"pair{i} = {color}" for i, color in enumerate(PAIR_COLORS.values(), 1))
    chroms = ";".join(COMMON_ORDER)
    config = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
{colors}
species_pzi = 0,114,178
species_psp = 0,158,115
species_ppa = 213,94,0
species_pcam = 117,99,177
species_ppe = 166,84,118
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {chroms}
<image>
background = white
dir = figures
file = comparative_pzi_psp_ppa_circos.png
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
radius = 0.800r
thickness = 7.06p
fill = yes
stroke_color = white
stroke_thickness = 1p
show_label = yes
label_radius = 1.045r
label_size = 22p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 0.800r
color = grey
thickness = 1p
multiplier = 1e-6
format = %d
<tick>
spacing = 50u
size = 3p
show_label = yes
label_size = 10p
label_offset = 2p
</tick>
</ticks>
<plots>
{chr(10).join(plots)}
</plots>
<links>
{links}
</links>
"""
    (out_root / "conf").mkdir(parents=True, exist_ok=True)
    (out_root / "conf" / "circos_integrated.conf").write_text(config, encoding="ascii")


def write_manifest(out_root: Path, common_lengths: dict[str, int]) -> None:
    with (out_root / "data" / "integrated_manifest.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("item", "value"))
        writer.writerow(("coordinate_system", "common canonical C1A-C8B; coordinates normalized per chromosome"))
        writer.writerow(("species_layers", "PZI, PSP, PPA"))
        writer.writerow(("tracks_per_species", "gene_density, te_coverage, ncrna_density, gc_fraction"))
        writer.writerow(("link_layers", "PZI A-B synteny; PZI-PSP, PZI-PPA, PZI-PCAM and PZI-PPE JCVI reference projections"))
        writer.writerow(("external_link_projections", "data/external_reference_mapping.tsv"))
        for chrom in COMMON_ORDER:
            writer.writerow((chrom, common_lengths[chrom]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--circos-root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--external-anchor-root", type=Path, default=EXTERNAL_ANCHOR_ROOT)
    parser.add_argument("--external-bed-root", type=Path, default=EXTERNAL_BED_ROOT)
    parser.add_argument("--work-anchor-root", type=Path, default=WORK_ANCHOR_ROOT)
    parser.add_argument("--work-bed-root", type=Path, default=WORK_BED_ROOT)
    args = parser.parse_args()
    circos_root = args.circos_root.resolve()
    out_root = args.out_root.resolve()
    data = out_root / "data"
    data.mkdir(parents=True, exist_ok=True)
    lengths = all_lengths(circos_root)
    remap_tracks(circos_root, data, lengths)
    remap_links(circos_root, data, lengths)
    remap_external_links(
        data,
        lengths,
        circos_root.resolve(),
        args.external_anchor_root.resolve(),
        args.external_bed_root.resolve(),
        args.work_anchor_root.resolve(),
        args.work_bed_root.resolve(),
    )
    write_karyotype(data, lengths)
    write_config(out_root, lengths)
    write_manifest(out_root, lengths)
    print(f"[OK] integrated common axes: {len(COMMON_ORDER)} chromosomes")
    print(f"[OK] output: {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
