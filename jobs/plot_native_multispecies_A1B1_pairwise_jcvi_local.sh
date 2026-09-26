#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare and render a PCAM-ordered, adjacent-track native JCVI ribbon figure.
# Run in WSL from the repository root.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR="${ROOT}/downstream_analysis/local_A1B1_ribbon_inputs"
ORDER_TABLE="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/multispecies_A1B1_pcam_aligned/chromosome_order_orientation.tsv"
OUT_DIR="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_multispecies_A1B1_pcam_pairwise"
PREPARE="${ROOT}/scripts/prepare_native_multispecies_A1B1_pairwise_jcvi.py"
RUNTIME_SHIM="${ROOT}/tmp/jcvi_native_runtime"
JCVI_PYTHON="${JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python3}"

require_file() { [[ -s "$1" ]] || { echo "[ERROR] Missing or empty file: $1" >&2; exit 1; }; }
require_file "${PREPARE}"
require_file "${ORDER_TABLE}"
require_file "${JCVI_PYTHON}"

"${JCVI_PYTHON}" "${PREPARE}" \
    --input-dir "${INPUT_DIR}" \
    --order-table "${ORDER_TABLE}" \
    --outdir "${OUT_DIR}" \
    --min-score 0 \
    --max-blocks-per-homolog 0

cd "${OUT_DIR}"
for format in png pdf svg; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        multispecies_A1B1_pairwise.seqids multispecies_A1B1_pairwise.layout \
        --figsize=18x24 --format="${format}" --notex --keep-chrlabels \
        --chrstyle=roundrect --shadestyle=curve \
        --outfile="multispecies_A1B1_pcam.pairwise.native.${format}"
done

for format in png pdf svg; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        multispecies_A1B1_pairwise.seqids multispecies_A1B1_pairwise_basepair.layout \
        --figsize=18x24 --format="${format}" --notex --nocircles --basepair \
        --chrstyle=roundrect --shadestyle=curve \
        --outfile="multispecies_A1B1_pcam.pairwise.native.basepair.${format}"
done

echo "[OK] Native JCVI pairwise figure: ${OUT_DIR}"
echo "[OK] Base-pair-scaled figure: ${OUT_DIR}/multispecies_A1B1_pcam.pairwise.native.basepair.png"
