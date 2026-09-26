#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT_ROOT}/0903}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/haploid_repartition/circos/pzi}"
GFF3="${GFF3:-${ROOT_DIR}/annotation/ncrna_clean/pzi/final/pzi.ncRNA.final.gff3}"
MAP="${MAP:-${ROOT_DIR}/haploid_repartition/maps/pzi_repartition_map.tsv}"
PREPARE="${PREPARE:-${PROJECT_ROOT}/scripts/prepare_pzi_ncrna_circos.py}"
LEGEND="${LEGEND:-${PROJECT_ROOT}/scripts/add_0903_circos_track_legend.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
CIRCOS_ENV="${CIRCOS_ENV:-circos}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

req "$GFF3" "cleaned PZI ncRNA GFF3"
req "$MAP" "PZI repartition map"
req "$PREPARE" "PZI ncRNA CircOS preparation script"
req "$LEGEND" "CircOS legend script"
req "$CONDA_SH" "conda initialization script"
req "$OUT_DIR/data/karyotype.txt" "existing PZI CircOS karyotype"
req "$OUT_DIR/conf/circos.conf" "existing PZI CircOS configuration"

mkdir -p "$OUT_DIR/data" "$OUT_DIR/conf" "$OUT_DIR/figures" "${ROOT_DIR}/logs"
source "$CONDA_SH"
conda activate "$CIRCOS_ENV"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${CIRCOS_ENV}"
command -v circos >/dev/null 2>&1 || die "circos is unavailable in ${CIRCOS_ENV}"

python "$PREPARE" \
    --gff3 "$GFF3" \
    --map "$MAP" \
    --karyotype "$OUT_DIR/data/karyotype.txt" \
    --out-data "$OUT_DIR/data" \
    --config "$OUT_DIR/conf/circos.conf" \
    --output-config "$OUT_DIR/conf/circos_ncrna.conf" \
    --output-stem pzi_ab_circos_ncrna

cd "$OUT_DIR"
rm -f figures/pzi_ab_circos_ncrna.png figures/pzi_ab_circos_ncrna.svg figures/pzi_ab_circos_ncrna.pdf
circos -conf conf/circos_ncrna.conf > "${OUT_DIR}/circos_ncrna_render.log" 2>&1
req "$OUT_DIR/figures/pzi_ab_circos_ncrna.png" "PZI ncRNA CircOS PNG"
req "$OUT_DIR/figures/pzi_ab_circos_ncrna.svg" "PZI ncRNA CircOS SVG"
grep -q '</svg>' "$OUT_DIR/figures/pzi_ab_circos_ncrna.svg" || die "PZI ncRNA SVG is incomplete"
python "$LEGEND" "$OUT_DIR/figures/pzi_ab_circos_ncrna.svg" \
    --png "$OUT_DIR/figures/pzi_ab_circos_ncrna.png" --include-ncrna

if command -v rsvg-convert >/dev/null 2>&1; then
    rsvg-convert -f pdf -o "$OUT_DIR/figures/pzi_ab_circos_ncrna.pdf" \
        "$OUT_DIR/figures/pzi_ab_circos_ncrna.svg"
elif command -v inkscape >/dev/null 2>&1; then
    inkscape "$OUT_DIR/figures/pzi_ab_circos_ncrna.svg" \
        --export-type=pdf --export-filename="$OUT_DIR/figures/pzi_ab_circos_ncrna.pdf"
else
    die "Need rsvg-convert or inkscape to convert the final SVG to PDF"
fi
req "$OUT_DIR/figures/pzi_ab_circos_ncrna.pdf" "PZI ncRNA CircOS PDF"

echo "[OK] PZI ncRNA CircOS outputs: $OUT_DIR/figures/pzi_ab_circos_ncrna.{png,svg,pdf}"
