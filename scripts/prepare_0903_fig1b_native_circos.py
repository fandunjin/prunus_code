#!/usr/bin/env python3
"""Prepare a native Circos Figure-1b style multispecies synteny plot.

This writes only Circos inputs.  The alignments are the existing native JCVI
``blocks.simple`` files under ``0903/native_multispecies_A1B1_pcam_pairwise``;
no new alignment is inferred here.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


AXIS = tuple([f"C{i}A" for i in range(1, 9)] + [f"C{i}B" for i in range(1, 9)])
PCAM_GROUPS = ("C1", "C3", "C7", "C8", "C4", "C5", "C6", "C2")
PAIR = {f"C{i}": f"pair{i}" for i in range(1, 9)}
EDGES = (("pzi", "psp"), ("psp", "ppa"), ("ppa", "pcam"), ("pcam", "ppe"), ("ppe", "gtr"), ("gtr", "rch"))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("0903"))
    p.add_argument("--native", type=Path, default=Path("0903/native_multispecies_A1B1_pcam_pairwise"))
    p.add_argument("--pzi-circos", type=Path, default=Path("0903/haploid_repartition/circos/pzi"))
    p.add_argument("--out", type=Path, default=Path("0903/haploid_repartition/circos/fig1b_native"))
    return p.parse_args()


def canonical(value: str) -> str | None:
    m = re.search(r"C0*(\d+)_?([AB])", value)
    return f"C{int(m.group(1))}{m.group(2)}" if m else None


def read_bed(path: Path):
    genes = {}
    lengths = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as h:
        for raw in h:
            f = raw.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            genes.setdefault(f[3], (f[0], s, e)); lengths[f[0]] = max(lengths[f[0]], e)
    return genes, dict(lengths)


def parse_karyotype(path: Path):
    out = {}
    if not path.is_file():
        return out
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        f = raw.split()
        if len(f) >= 6 and f[0] == "chr":
            out[f[2]] = int(f[5])
    return out


def parse_blocks(path: Path):
    with path.open(encoding="utf-8", errors="replace") as h:
        for raw in h:
            line = raw.strip()
            if not line or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            f = line.split()
            if len(f) == 6:
                try:
                    yield f[0], f[1], f[2], f[3], int(float(f[4])), f[5]
                except ValueError:
                    pass


def seq_number(value: str) -> int:
    m = re.search(r"\d+", value)
    return int(m.group()) if m else 99


def build_maps(native: Path, beds: dict[str, tuple[dict, dict]]):
    maps = {}
    for species in ("pzi", "psp", "ppa"):
        maps[species] = {seq: canonical(seq) for seq in beds[species][1] if canonical(seq) in AXIS}
    for species in ("pcam", "ppe"):
        seqids = sorted(beds[species][1], key=seq_number)
        maps[species] = {seq: f"{PCAM_GROUPS[i]}B" for i, seq in enumerate(seqids[:8])}

    # GTR chromosomes are assigned from the ppe--gtr blocks.  The ppe track is
    # on B, but GTR is shown on the A side of the common PZI axis.
    ppe_genes, _ = beds["ppe"]; gtr_genes, _ = beds["gtr"]
    votes = defaultdict(lambda: defaultdict(int))
    for a, b, c, d, score, _ in parse_blocks(native / "anchors" / "ppe.gtr.blocks.simple"):
        if a not in ppe_genes or c not in gtr_genes:
            continue
        pseq, gseq = ppe_genes[a][0], gtr_genes[c][0]
        pc = maps["ppe"].get(pseq)
        if pc:
            votes[gseq][pc[:2] + "A"] += max(1, score)
    gtr_order = sorted(beds["gtr"][1], key=seq_number)
    maps["gtr"] = {seq: max(votes[seq], key=votes[seq].get) if votes[seq] else f"{PCAM_GROUPS[i % 8]}A" for i, seq in enumerate(gtr_order)}

    # RCH is placed on the A side using the GTR--RCH block votes.
    rch_genes, _ = beds["rch"]; votes_r = defaultdict(lambda: defaultdict(int))
    for a, b, c, d, score, _ in parse_blocks(native / "anchors" / "gtr.rch.blocks.simple"):
        if a not in gtr_genes or c not in rch_genes:
            continue
        gc = maps["gtr"].get(gtr_genes[a][0])
        if gc:
            votes_r[rch_genes[c][0]][gc] += max(1, score)
    maps["rch"] = {seq: max(votes_r[seq], key=votes_r[seq].get) if votes_r[seq] else f"{PCAM_GROUPS[i % 8]}A" for i, seq in enumerate(sorted(beds["rch"][1], key=seq_number))}
    return maps


def source_lengths(pzi: Path, species: str, beds):
    if species in ("pzi", "psp", "ppa"):
        path = pzi.parent / species / "data" / "karyotype.txt"
        rows = parse_karyotype(path)
        result = {}
        for seq, length in rows.items():
            c = canonical(seq)
            if not c:
                continue
            if species == "ppa" and seq.startswith(("A_", "B_")):
                result[c] = length
            elif species != "ppa" and seq.startswith(("A1_", "B1_")):
                result[c] = length
        return result
    return beds[species][1]


def remap(x: int, old: int, new: int) -> int:
    return max(0, min(new, round(x * new / max(1, old))))


def write_axis(data: Path, lengths: dict[str, int]):
    with (data / "karyotype.txt").open("w", encoding="ascii") as h:
        for c in AXIS:
            h.write(f"chr - {c} {c} 0 {lengths[c]} {PAIR[c[:2]]}\n")


def write_track(data: Path, pzi: Path, name: str):
    src = pzi / "data" / f"{name}.txt"; dst = data / f"pzi_{name}.txt"; rows = []
    if src.is_file():
        for raw in src.read_text(encoding="utf-8", errors="replace").splitlines():
            f = raw.split()
            if len(f) < 4:
                continue
            c = canonical(f[0].split("_", 1)[1]) if "_" in f[0] else canonical(f[0])
            if c and (f[0].startswith("A1_") or f[0].startswith("B1_")):
                rows.append(f"{c} {f[1]} {f[2]} {f[3]}")
    dst.write_text("\n".join(rows) + "\n", encoding="ascii")


def write_bands(data: Path, maps, beds, lengths):
    for species, mapping in maps.items():
        rows = []
        for seq, c in mapping.items():
            if not c or c not in lengths or seq not in beds[species][1]:
                continue
            src = beds[species][1][seq]; axis_len = lengths[c]
            width = min(axis_len * 0.94, src); center = axis_len / 2
            rows.append(f"{c} {round(center-width/2)} {round(center+width/2)} color=black")
        (data / f"{species}_bands.txt").write_text("\n".join(rows) + "\n", encoding="ascii")


def write_links(data: Path, native: Path, beds, maps, lengths):
    counts = {}
    for upper, lower in EDGES:
        filename = f"{upper}.{lower}.blocks.simple"; rows = []
        ug, ul = beds[upper]; lg, ll = beds[lower]
        for a, b, c, d, score, orient in parse_blocks(native / "anchors" / filename):
            if not all(x in ug for x in (a, b)) or not all(x in lg for x in (c, d)):
                continue
            q, q2, s, s2 = ug[a], ug[b], lg[c], lg[d]
            if q[0] != q2[0] or s[0] != s2[0]:
                continue
            qc, sc = maps[upper].get(q[0]), maps[lower].get(s[0])
            if not qc or not sc or qc not in lengths or sc not in lengths:
                continue
            rows.append(f"{qc} {remap(min(q[1],q2[1]),ul[q[0]],lengths[qc])} {remap(max(q[2],q2[2]),ul[q[0]],lengths[qc])} {sc} {remap(min(s[1],s2[1]),ll[s[0]],lengths[sc])} {remap(max(s[2],s2[2]),ll[s[0]],lengths[sc])} color={PAIR.get(qc[:2],'pair8')}_a70,z=20,twist={'yes' if orient == '-' else 'no'}")
        (data / f"links_{upper}_{lower}.txt").write_text("\n".join(rows) + "\n", encoding="ascii"); counts[f"{upper}-{lower}"] = len(rows)
    return counts


def write_internal(data: Path, pzi: Path, lengths):
    rows = []
    for raw in (pzi / "data" / "links_ab.txt").read_text(encoding="utf-8", errors="replace").splitlines():
        f = raw.split()
        if len(f) < 6 or not f[0].startswith("A1_") or not f[3].startswith("B1_"):
            continue
        q, s = canonical(f[0].split("_",1)[1]), canonical(f[3].split("_",1)[1])
        if q and s:
            rows.append(f"{q} {f[1]} {f[2]} {s} {f[4]} {f[5]} color={PAIR[q[:2]]}_a38,z=100,twist=no")
    (data / "links_pzi_A1_B1.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_config(out: Path, maxima):
    conf = f'''<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
pair1 = 0,114,178
pair2 = 230,159,0
pair3 = 0,158,115
pair4 = 213,94,0
pair5 = 204,121,167
pair6 = 86,180,233
pair7 = 190,168,0
pair8 = 80,80,80
gene_track = 85,90,95
te_track = 112,135,52
gypsy_track = 196,91,82
copia_track = 76,120,168
tir_track = 126,88,160
helitron_track = 40,150,145
gc_track = 42,88,108
ncrna_track = 111,158,158
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {";".join(AXIS)}
<image>
background = white
dir = figures
file = fig1b_integrated_pzi_multispecies.png
png = yes
svg = yes
radius = 900p
angle_offset = 90
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.004r
</spacing>
radius = 0.620r
thickness = 5p
fill = yes
stroke_color = black
stroke_thickness = 0.5p
show_label = yes
label_radius = 1.030r
label_size = 8p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 0.620r
color = grey
thickness = 0.5p
multiplier = 1e-6
format = %d
<tick>
spacing = 50u
size = 4p
show_label = yes
label_size = 6p
label_offset = 2p
</tick>
</ticks>
<plots>
<plot>
type = heatmap
file = data/pzi_gene_density.txt
r0 = 0.575r
r1 = 0.610r
min = 0
max = {maxima.get('gene_density',1):.0f}
color = spectral-11-div
stroke_thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_te_coverage.txt
r0 = 0.525r
r1 = 0.560r
min = 0
max = 1
fill_color = te_track
color = te_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_gypsy_coverage.txt
r0 = 0.475r
r1 = 0.510r
min = 0
max = 0.32
fill_color = gypsy_track
color = gypsy_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_copia_coverage.txt
r0 = 0.425r
r1 = 0.460r
min = 0
max = 0.16
fill_color = copia_track
color = copia_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_tir_coverage.txt
r0 = 0.375r
r1 = 0.410r
min = 0
max = 0.41
fill_color = tir_track
color = tir_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_helitron_coverage.txt
r0 = 0.325r
r1 = 0.360r
min = 0
max = 0.06
fill_color = helitron_track
color = helitron_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/pzi_ncrna_density.txt
r0 = 0.275r
r1 = 0.310r
min = 0
max = 1000
fill_color = ncrna_track
color = ncrna_track
thickness = 0
</plot>
<plot>
type = line
file = data/pzi_gc_fraction.txt
r0 = 0.235r
r1 = 0.265r
min = 0.35
max = 0.56
color = gc_track
thickness = 1p
</plot>
</plots>
<highlights>
<highlight>
file = data/pzi_bands.txt
r0 = 580p
r1 = 588p
</highlight>
<highlight>
file = data/psp_bands.txt
r0 = 620p
r1 = 628p
</highlight>
<highlight>
file = data/ppa_bands.txt
r0 = 660p
r1 = 668p
</highlight>
<highlight>
file = data/pcam_bands.txt
r0 = 700p
r1 = 708p
</highlight>
<highlight>
file = data/ppe_bands.txt
r0 = 740p
r1 = 748p
</highlight>
<highlight>
file = data/gtr_bands.txt
r0 = 780p
r1 = 788p
</highlight>
<highlight>
file = data/rch_bands.txt
r0 = 820p
r1 = 828p
</highlight>
</highlights>
<links>
<link>
file = data/links_pzi_A1_B1.txt
radius = 0.220r
bezier_radius = 0.070r
ribbon = yes
thickness = 0.6p
crest = 0.45
</link>
<link>
file = data/links_pzi_psp.txt
radius1 = 588p
radius2 = 620p
bezier_radius = 604p
bezier_radius_purity = 1
flat = yes
ribbon = yes
thickness = 0.6p
crest = 1
</link>
<link>
file = data/links_psp_ppa.txt
radius1 = 628p
radius2 = 660p
bezier_radius = 644p
bezier_radius_purity = 1
flat = yes
ribbon = yes
thickness = 0.6p
crest = 1
</link>
<link>
file = data/links_ppa_pcam.txt
radius1 = 668p
radius2 = 700p
bezier_radius = 684p
bezier_radius_purity = 1
flat = yes
ribbon = no
thickness = 1.2p
crest = 1
</link>
<link>
file = data/links_pcam_ppe.txt
radius1 = 708p
radius2 = 740p
bezier_radius = 724p
bezier_radius_purity = 1
flat = yes
ribbon = yes
thickness = 0.6p
crest = 1
</link>
<link>
file = data/links_ppe_gtr.txt
radius1 = 748p
radius2 = 780p
bezier_radius = 764p
bezier_radius_purity = 1
flat = yes
ribbon = no
thickness = 1.2p
crest = 1
</link>
<link>
file = data/links_gtr_rch.txt
radius1 = 788p
radius2 = 820p
bezier_radius = 804p
bezier_radius_purity = 1
flat = yes
ribbon = yes
thickness = 0.6p
crest = 1
</link>
</links>
'''
    (out / "conf").mkdir(parents=True, exist_ok=True); (out / "conf" / "circos.conf").write_text(conf, encoding="ascii")


def main():
    a = parse_args(); native, pzi, out = a.native.resolve(), a.pzi_circos.resolve(), a.out.resolve(); data = out / "data"; data.mkdir(parents=True, exist_ok=True); (out / "figures").mkdir(parents=True, exist_ok=True)
    beds = {s: read_bed(native / "bed" / f"{s}.bed") for s in ("pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch")}
    lengths = source_lengths(pzi, "pzi", beds)
    for c in AXIS: lengths.setdefault(c, max([beds["pzi"][1].get(c, 1)] + [1]))
    maps = build_maps(native, beds)
    write_axis(data, lengths)
    maxima = {}
    for name in ("gene_density", "te_coverage", "gypsy_coverage", "copia_coverage", "tir_coverage", "helitron_coverage", "gc_fraction", "ncrna_density"):
        write_track(data, pzi, name)
        vals = []
        for raw in (data / f"pzi_{name}.txt").read_text().splitlines():
            f = raw.split()
            if len(f) >= 4:
                try: vals.append(float(f[3]))
                except ValueError: pass
        maxima[name] = max(vals or [1])
    write_bands(data, maps, beds, lengths)
    internal = write_internal(data, pzi, lengths); counts = write_links(data, native, beds, maps, lengths); write_config(out, maxima)
    with (out / "figure_manifest.tsv").open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(("item", "value")); w.writerow(("canonical_axis", ",".join(AXIS))); w.writerow(("pzi_A1_B1_links", internal)); w.writerow(("external_link_geometry", "same-side=ribbon; cross-side=PPA-PCAM,PPE-GTR=thin annular curve"))
        for k, v in counts.items(): w.writerow((k, v))
        for s, mp in maps.items():
            for seq, c in mp.items(): w.writerow((f"{s}:{seq}", c))
    (out / "README_methods.txt").write_text("Native Circos Figure-1b style plot. PZI A1/B1 are the inner axes with existing PZI density, TE, GC and ncRNA tracks. Outer black chromosome bands are centered on the common C1A-C8B axes; band widths use source chromosome bp lengths. The seven species bands occupy explicit outer radii 580-828p. Same-side adjacent comparisons (PZI-PSP, PSP-PPA, PCAM-PPE and GTR-RCH) are filled ribbons constrained to their adjacent radial intervals. Cross-side comparisons (PPA-PCAM and PPE-GTR) retain all native JCVI links but are rendered as thin curves in their adjacent radial intervals to avoid a central filled chord.\n", encoding="ascii")
    print(f"[OK] prepared {out}"); print(f"[OK] PZI internal links: {internal}"); [print(f"[OK] {k}: {v} links") for k, v in counts.items()]


if __name__ == "__main__":
    main()
