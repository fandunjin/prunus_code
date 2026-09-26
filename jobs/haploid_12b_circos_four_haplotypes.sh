#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/circos}"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_pzi_four_haplotype_circos_jcvi.py"
GUIDE_SCRIPT="${ROOT_DIR}/scripts/add_circos_gap_track_guide.py"
QUARTET_TSV="${ROOT_DIR}/haploid/final_haplotypes/details/final_quartet_assignment.tsv"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
CIRCOS_ENV="${CIRCOS_ENV:-circos}"

WINDOW="${WINDOW:-1000000}"
FORCE_PREP="${FORCE_PREP:-0}"
FORCE_PLOT="${FORCE_PLOT:-0}"

JCVI_DIR="${OUT_DIR}/work/jcvi"
CURRENT_SIGNATURE="${OUT_DIR}/.circos_inputs.current"
SAVED_SIGNATURE="${OUT_DIR}/.circos_inputs.ok"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

valid_svg() {
    [[ -s "$1" ]] && grep -q '</svg>' "$1"
}

mkdir -p "${OUT_DIR}/data" "${OUT_DIR}/conf" "${OUT_DIR}/figures" "${ROOT_DIR}/logs"
require_file "${CONDA_SH}"
require_file "${PREP_SCRIPT}"
require_file "${GUIDE_SCRIPT}"
require_file "${QUARTET_TSV}"
require_file "${HAP_RENAME_DIR}/chromosome_rename_map.tsv"
require_file "${HAP_RENAME_DIR}/checksums.md5"
for hap in A1 A2 B1 B2; do
    require_file "${HAP_RENAME_DIR}/genome/pzi_hap${hap}.fa"
    require_file "${HAP_RENAME_DIR}/gff/pzi_hap${hap}.gff3"
    require_file "${HAP_RENAME_DIR}/cds/pzi_hap${hap}.cds.fa"
    require_file "${HAP_RENAME_DIR}/pep/pzi_hap${hap}.pep.fa"
done

for hap in A1 A2 B1 B2; do
    require_file "${OUT_DIR}/work/jcvi_inputs/hap${hap}.bed"
done
for pair in A1_A2 B1_B2 A1_B1 A2_B2; do
    require_file "${JCVI_DIR}/${pair}/blocks.simple"
done

{
    printf 'window=%s\n' "${WINDOW}"
    cksum \
        "${PREP_SCRIPT}" \
        "${GUIDE_SCRIPT}" \
        "${QUARTET_TSV}" \
        "${HAP_RENAME_DIR}/chromosome_rename_map.tsv" \
        "${HAP_RENAME_DIR}/checksums.md5" \
        "${JCVI_DIR}/A1_A2/blocks.simple" \
        "${JCVI_DIR}/B1_B2/blocks.simple" \
        "${JCVI_DIR}/A1_B1/blocks.simple" \
        "${JCVI_DIR}/A2_B2/blocks.simple"
} > "${CURRENT_SIGNATURE}"

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${CIRCOS_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${CIRCOS_ENV}"
command -v circos >/dev/null 2>&1 || die "circos is unavailable in ${CIRCOS_ENV}"

echo "=== P. zippeliana four-haplotype Circos rendering ==="
echo "Root       : ${ROOT_DIR}"
echo "Hap package: ${HAP_RENAME_DIR}"
echo "Output     : ${OUT_DIR}"
echo "JCVI blocks: ${JCVI_DIR}"
echo "Circos env : ${CIRCOS_ENV}"
echo "Window     : ${WINDOW} bp"

PREP_RAN=0
if [[ "${FORCE_PREP}" == "1" ]] || ! cmp -s "${CURRENT_SIGNATURE}" "${SAVED_SIGNATURE}"; then
    python "${PREP_SCRIPT}" \
        --root "${ROOT_DIR}" \
        --hap-rename-dir "${HAP_RENAME_DIR}" \
        --out-dir "${OUT_DIR}" \
        --jcvi-dir "${JCVI_DIR}" \
        --window "${WINDOW}"
    cp -f "${CURRENT_SIGNATURE}" "${SAVED_SIGNATURE}"
    PREP_RAN=1
else
    echo "[SKIP] Existing Circos inputs match the JCVI blocks and window size"
fi

require_file "${OUT_DIR}/conf/circos.conf"
require_file "${OUT_DIR}/data/links_background.txt"
require_file "${OUT_DIR}/data/links_highlight.txt"
require_file "${OUT_DIR}/data/gypsy_coverage.txt"
require_file "${OUT_DIR}/data/copia_coverage.txt"

PNG="${OUT_DIR}/figures/pzi_four_haplotypes_jcvi_circos.png"
SVG="${OUT_DIR}/figures/pzi_four_haplotypes_jcvi_circos.svg"
PDF="${OUT_DIR}/figures/pzi_four_haplotypes_jcvi_circos.pdf"

PLOT_RAN=0
NEED_PLOT=0
if [[ "${FORCE_PLOT}" == "1" || "${PREP_RAN}" == "1" || ! -s "${PNG}" ]]; then
    NEED_PLOT=1
elif ! valid_svg "${SVG}"; then
    NEED_PLOT=1
fi

if [[ "${NEED_PLOT}" == "1" ]]; then
    # A previous Circos failure can leave a non-empty but truncated SVG.
    rm -f "${PNG}" "${SVG}"
    cd "${OUT_DIR}"
    circos -conf conf/circos.conf
    PLOT_RAN=1
else
    echo "[SKIP] Existing Circos PNG and SVG"
fi
require_file "${PNG}"
require_file "${SVG}"
valid_svg "${SVG}" || die "Circos produced an incomplete SVG: ${SVG}"
require_file "${OUT_DIR}/data/track_layout.tsv"

# Circos cannot natively place a free radial guide in a blank angular sector.
# Add the A-E guide after the SVG is rendered; this is idempotent on resumes.
python "${GUIDE_SCRIPT}" \
    --svg "${SVG}" \
    --layout "${OUT_DIR}/data/track_layout.tsv" \
    --ideogram-radius-factor 0.940 \
    --ideogram-thickness-p 22
require_file "${SVG}"

if command -v rsvg-convert >/dev/null 2>&1; then
    # Keep the raster preview consistent with the labelled SVG/PDF.
    rsvg-convert -f png -o "${PNG}" "${SVG}"
    require_file "${PNG}"
    if [[ "${FORCE_PLOT}" == "1" || "${PREP_RAN}" == "1" || "${PLOT_RAN}" == "1" || ! -s "${PDF}" ]]; then
        rsvg-convert -f pdf -o "${PDF}" "${SVG}"
    fi
else
    echo "[WARN] rsvg-convert is unavailable; PDF conversion was skipped."
fi

echo
echo "=== Circos outputs ==="
ls -lh "${OUT_DIR}/figures"/pzi_four_haplotypes_jcvi_circos.*
echo "Track table : ${OUT_DIR}/data/window_tracks.tsv"
echo "JCVI summary: ${OUT_DIR}/data/jcvi_block_summary.tsv"
echo "Link audit  : ${OUT_DIR}/data/link_coverage_by_chromosome.tsv"
echo "Layout scale: ${OUT_DIR}/data/layout_scaling.tsv"
echo "Chr aliases : ${OUT_DIR}/data/chromosome_rename_map.tsv"
echo "Orientation : ${OUT_DIR}/data/display_orientation.tsv"
echo "Copy1 pairs : ${OUT_DIR}/data/copy1_pair_map.tsv"
echo "Track layout: ${OUT_DIR}/data/track_layout.tsv"
