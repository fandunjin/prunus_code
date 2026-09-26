#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

# Draw a native JCVI ribbon with one PZI copy1 chromosome per horizontal track.
# The track order is controlled by ORDER and is preserved exactly.

ROOT_DIR="${ROOT_DIR:-$(pwd)}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_selected_vertical_ribbon_jcvi_local}"
ORDER="${ORDER:-1 13 2 15 5 11 16 6 3}"
JCVI_PYTHON="${JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python}"
RUNTIME_SHIM="${RUNTIME_SHIM:-${ROOT_DIR}/tmp/jcvi_native_runtime}"
SOURCE_BED="${SOURCE_BED:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/pzi_copy1.bed}"
SOURCE_ANCHORS="${SOURCE_ANCHORS:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/pzi_copy1.pzi_copy1.anchors.simple}"
FIGSIZE="${FIGSIZE:-10x20}"

die() { echo "[ERROR] $*" >&2; exit 1; }
req() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

req "$JCVI_PYTHON"
req "$SOURCE_BED"
req "$SOURCE_ANCHORS"
[[ -d "$RUNTIME_SHIM" ]] || die "Missing JCVI runtime shim: $RUNTIME_SHIM"

BED_DIR="$OUT_DIR/bed"
ANCHOR_DIR="$OUT_DIR/anchors"
FIG_DIR="$OUT_DIR/figures"
mkdir -p "$BED_DIR" "$ANCHOR_DIR" "$FIG_DIR"

export OUT_DIR ORDER SOURCE_BED SOURCE_ANCHORS BED_DIR ANCHOR_DIR FIG_DIR
"$JCVI_PYTHON" - <<'PY'
import os
from pathlib import Path

order = [int(x) for x in os.environ["ORDER"].split()]
chroms = [f"chr{x:02d}" for x in order]
source_bed = Path(os.environ["SOURCE_BED"])
bed_dir = Path(os.environ["BED_DIR"])
anchor_dir = Path(os.environ["ANCHOR_DIR"])
source_anchors = Path(os.environ["SOURCE_ANCHORS"])

gene_chr = {}
by_chr = {chrom: [] for chrom in chroms}
for raw in source_bed.read_text(encoding="utf-8", errors="replace").splitlines():
    if not raw.strip() or raw.startswith("#"):
        continue
    fields = raw.split("\t")
    if len(fields) < 6:
        continue
    chrom, gene = fields[0], fields[3]
    gene_chr[gene] = chrom
    if chrom in by_chr:
        by_chr[chrom].append(raw)

for chrom in chroms:
    rows = by_chr[chrom]
    if not rows:
        raise SystemExit(f"No BED genes found for {chrom}")
    (bed_dir / f"{chrom}.bed").write_text("\n".join(rows) + "\n", encoding="utf-8")

pair_rows = {(chroms[i], chroms[i + 1]): [] for i in range(len(chroms) - 1)}
for raw in source_anchors.read_text(encoding="utf-8", errors="replace").splitlines():
    line = raw.strip()
    if not line or line.startswith("#") or line.startswith("StartGeneA"):
        continue
    if "*" in line:
        line = line.split("*", 1)[1].strip()
    fields = line.split()
    if len(fields) != 6:
        continue
    a, b, c, d, score, orient = fields
    qchr, schr = gene_chr.get(a), gene_chr.get(c)
    if not qchr or not schr or qchr == schr:
        continue
    key = (qchr, schr)
    reverse_key = (schr, qchr)
    if key in pair_rows:
        pair_rows[key].append(fields)
    elif reverse_key in pair_rows:
        pair_rows[reverse_key].append([c, d, a, b, score, orient])

edge_pairs = []
for (left, right), rows in pair_rows.items():
    path = anchor_dir / f"{left}.{right}.blocks.simple"
    if not rows:
        continue
    path.write_text("\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")
    edge_pairs.append((left, right))

seqids = Path(os.environ["OUT_DIR"]) / "selected_vertical.seqids"
seqids.write_text("\n".join(chroms) + "\n", encoding="utf-8")

layout = Path(os.environ["OUT_DIR"]) / "selected_vertical.layout"
lines = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
n = len(chroms)
for i, chrom in enumerate(chroms):
    y = 0.96 - i * (0.88 / max(1, n - 1))
    lines.append(f" {y:.5f}, .10, .90, 0, #4C5660, {chrom}, center, bed/{chrom}.bed, center")
lines.append("# edges")
for i in range(n - 1):
    left, right = chroms[i], chroms[i + 1]
    if (left, right) in edge_pairs:
        lines.append(f" e, {i}, {i+1}, anchors/{left}.{right}.blocks.simple")
layout.write_text("\n".join(lines) + "\n", encoding="utf-8")

basepair = Path(os.environ["OUT_DIR"]) / "selected_vertical_basepair.layout"
basepair.write_text(layout.read_text(encoding="utf-8"), encoding="utf-8")

summary = Path(os.environ["OUT_DIR"]) / "order.tsv"
summary.write_text("rank\tcopy1_chr\n" + "\n".join(f"{i+1}\t{chrom}" for i, chrom in enumerate(chroms)) + "\n", encoding="utf-8")
PY

cd "$OUT_DIR"
for format in png pdf svg; do
    PYTHONPATH="$RUNTIME_SHIM${PYTHONPATH:+:${PYTHONPATH}}" \
        "$JCVI_PYTHON" -m jcvi.graphics.karyotype \
        selected_vertical.seqids selected_vertical.layout \
        --figsize="$FIGSIZE" --format="$format" --notex --nocircles \
        --keep-chrlabels --chrstyle=roundrect --shadestyle=curve \
        --outfile="$FIG_DIR/pzi_copy1_selected_vertical_ribbon.${format}"
done

echo "[OK] PZI copy1 vertical ribbon written to: $FIG_DIR"
echo "[OK] Order: $ORDER"
