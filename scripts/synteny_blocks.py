from pathlib import Path
import csv, collections, json
import numpy as np

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project-root", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
O = args.out
R = args.project_root / "0903"
old = R / "native_multispecies_A1B1_pcam_pairwise"

colors = [
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#BEA800",
    "#505050",
]
beds = {}
for s in ["pzi", "psp", "ppa"]:
    beds[s] = {
        v[3]: (v[0].split("__")[-1], int(v[1]), int(v[2]))
        for v in [x.split() for x in (O / f"{s}.bed").read_text().splitlines()]
    }
allblocks = []
newpairs = {}
for left, right in [("pzi", "psp"), ("psp", "ppa"), ("pzi", "ppa")]:
    groups = []
    current = []
    for line in (O / f"{left}.{right}.anchors").open():
        if line.startswith("#"):
            if current:
                groups.append(current)
            current = []
        else:
            v = line.split()
            if len(v) >= 2:
                current.append((v[0], v[1]))
    if current:
        groups.append(current)
    candidates = collections.defaultdict(list)
    for i, g in enumerate(groups):
        a = [beds[left][x] for x, y in g]
        b = [beds[right][y] for x, y in g]
        if len({v[0] for v in a}) != 1 or len({v[0] for v in b}) != 1:
            continue
        lc, rc = a[0][0], b[0][0]
        # Match structural chromosome counterparts, as in the original interspecific display.
        if lc != rc:
            continue
        xs = np.array([v[1] for v in a])
        ys = np.array([v[1] for v in b])
        ori = "+" if np.cov(xs, ys)[0, 1] >= 0 else "-"
        q0 = min(a, key=lambda v: v[1])
        q1 = max(a, key=lambda v: v[2])
        s0 = min(b, key=lambda v: v[1])
        s1 = max(b, key=lambda v: v[2])
        rec = dict(
            pair=left + "-" + right,
            left_species=left,
            right_species=right,
            left_chromosome=lc,
            left_start0=q0[1],
            left_end=q1[2],
            right_chromosome=rc,
            right_start0=s0[1],
            right_end=s1[2],
            anchors=len(g),
            orientation=ori,
            source_block=i + 1,
            left_start_gene=min(g, key=lambda z: beds[left][z[0]][1])[0],
            left_end_gene=max(g, key=lambda z: beds[left][z[0]][2])[0],
            right_start_gene=min(g, key=lambda z: beds[right][z[1]][1])[1],
            right_end_gene=max(g, key=lambda z: beds[right][z[1]][2])[1],
            method="DIAMOND cscore>=0.70; JCVI scan min_size=4 dist=30; up to four largest blocks per matching chromosome",
        )
        candidates[lc].append(rec)
    keep = []
    for c, rs in candidates.items():
        # Largest anchor-supported blocks; reject almost completely redundant spans.
        selected = []
        for r in sorted(rs, key=lambda r: r["anchors"], reverse=True):
            if any(
                max(
                    0,
                    min(r["left_end"], s["left_end"])
                    - max(r["left_start0"], s["left_start0"]),
                )
                / (r["left_end"] - r["left_start0"])
                > 0.8
                and max(
                    0,
                    min(r["right_end"], s["right_end"])
                    - max(r["right_start0"], s["right_start0"]),
                )
                / (r["right_end"] - r["right_start0"])
                > 0.8
                for s in selected
            ):
                continue
            selected.append(r)
            if len(selected) == 4:
                break
        keep.extend(selected)
    newpairs[left, right] = keep
    allblocks += keep
    print(left, right, len(keep), "chromosomes", len(candidates))
# Add unchanged, existing non-focal adjacent comparisons with actual endpoint coordinates.
for s in ["ppa", "pcam", "ppe", "gtr", "rch"]:
    if s == "ppa":
        # Legacy blocks may name an alternative transcript, so retain full endpoint catalog for tracing.
        pass
    beds[s + "_legacy"] = {
        v[3]: (v[0], int(v[1]), int(v[2]))
        for v in [x.split() for x in (old / f"bed/{s}.bed").read_text().splitlines()]
    }
for left, right in [("ppa", "pcam"), ("pcam", "ppe"), ("ppe", "gtr"), ("gtr", "rch")]:
    rs = []
    for i, line in enumerate(
        (old / f"anchors/{left}.{right}.blocks.simple").read_text().splitlines()
    ):
        v = line.split("*")[-1].split()
        g1, g2, h1, h2, n, ori = v
        lb = beds[left + "_legacy"]
        rb = beds[right + "_legacy"]
        a, b = lb[g1], lb[g2]
        c, d = rb[h1], rb[h2]
        rs.append(
            dict(
                pair=left + "-" + right,
                left_species=left,
                right_species=right,
                left_chromosome=a[0],
                left_start0=min(a[1], b[1]),
                left_end=max(a[2], b[2]),
                right_chromosome=c[0],
                right_start0=min(c[1], d[1]),
                right_end=max(c[2], d[2]),
                anchors=int(n),
                orientation=ori,
                source_block=i + 1,
                left_start_gene=g1,
                left_end_gene=g2,
                right_start_gene=h1,
                right_end_gene=h2,
                method="Retained archived adjacent-species block",
            )
        )
    newpairs[left, right] = rs
    allblocks += rs
# Save separate all-data and displayed circos/riparian data.
with (O / "all_display_blocks.tsv").open("w") as f:
    wr = csv.DictWriter(f, fieldnames=list(allblocks[0]), delimiter="\t")
    wr.writeheader()
    wr.writerows(allblocks)
chromrows = list(csv.DictReader((O / "chromosomes.tsv").open(), delimiter="\t"))
lengths = collections.defaultdict(dict)
for r in chromrows:
    lengths[r["species"]][r["chromosome"]] = int(r["length_bp"])
for sp in ["pcam", "ppe", "gtr", "rch"]:
    for c, a, b in beds[sp + "_legacy"].values():
        lengths[sp][c] = max(lengths[sp].get(c, 0), b)
for r in allblocks:
    assert r["left_end"] <= lengths[r["left_species"]][r["left_chromosome"]], r
    assert r["right_end"] <= lengths[r["right_species"]][r["right_chromosome"]], r
# Circos links use current actual A1/B1 coordinates and direct pairwise blocks.
for pair, outname in [
    (("pzi", "psp"), "pzi_psp"),
    (("psp", "ppa"), "psp_ppa"),
    (("pzi", "ppa"), "ppa_pzi"),
]:
    rs = newpairs[pair]
    lines = []
    for i, r in enumerate(rs):
        # Circos ribbon twist follows the inferred block orientation.
        twist = "yes" if r["orientation"] == "-" else "no"
        lines.append(
            f"{r['left_species']}_{r['left_chromosome']} {r['left_start0']} {r['left_end']} {r['right_species']}_{r['right_chromosome']} {r['right_start0']} {r['right_end']} color=link_{outname},z={100+i},twist={twist}"
        )
    (O / f"circos/data/links_{outname}.txt").write_text("\n".join(lines) + "\n")
# Preserve original species and chromosome order; orient the new PZI homologs against PSP.
ss = ["pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"]
order = {}
flips = {}
for sp, line in zip(
    ss, (old / "multispecies_0903_A1B1_pcam.seqids").read_text().splitlines()
):
    order[sp] = [c.rstrip("-") for c in line.split(",")]
    flips[sp] = {c.rstrip("-"): c.endswith("-") for c in line.split(",")}
for c in order["pzi"]:
    r = max(
        [x for x in newpairs["pzi", "psp"] if x["left_chromosome"] == c],
        key=lambda z: z["anchors"],
    )
    flips["pzi"][c] = flips["psp"][c] ^ (r["orientation"] == "-")
for c in order["ppa"]:
    r = max(
        [x for x in newpairs["psp", "ppa"] if x["left_chromosome"] == c],
        key=lambda z: z["anchors"],
    )
    flips["ppa"][c] = flips["psp"][c] ^ (r["orientation"] == "-")
ss = ["pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"]
layout = []
for sp, y in zip(ss, [0.95, 0.81, 0.67, 0.53, 0.39, 0.25, 0.11]):
    for c in order[sp]:
        layout.append(
            dict(
                species=sp,
                chromosome=c,
                length_bp=lengths[sp][c],
                display_flipped=flips[sp][c],
                xstart=0,
                xend=0,
                y=y,
            )
        )
with (O / "riparian_layout.tsv").open("w") as f:
    wr = csv.DictWriter(f, fieldnames=list(layout[0]), delimiter="\t")
    wr.writeheader()
    wr.writerows(layout)
(O / "ribbon_validation.json").write_text(
    json.dumps(
        {
            "block_counts": {"-".join(k): len(v) for k, v in newpairs.items()},
            "all_endpoints_within_chromosome": True,
            "focal_inputs": "current representative loci, PZI/PSP A1+B1 and PPA A+B",
            "chromosome_count_by_species": {s: len(order[s]) for s in ss},
        },
        indent=2,
    )
)
