#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare and render a native JCVI karyotype figure with selected PZI self blocks.
# Run from WSL at the repository root.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREPARE="${ROOT}/scripts/prepare_native_pzi_self_highlight_jcvi.py"
OUT_DIR="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_pzi_self_highlight_ppa_psp_ppe_pcam_gtr"
RUNTIME_SHIM="${ROOT}/tmp/jcvi_native_runtime"
JCVI_PYTHON="${JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python3}"

[[ -s "${PREPARE}" ]] || { echo "[ERROR] Missing prepare script: ${PREPARE}" >&2; exit 1; }
[[ -x "${JCVI_PYTHON}" ]] || { echo "[ERROR] JCVI Python is not executable: ${JCVI_PYTHON}" >&2; exit 1; }

"${JCVI_PYTHON}" "${PREPARE}" \
    --outdir "${OUT_DIR}" \
    --min-score 0 \
    --max-blocks-per-homolog 0

cd "${OUT_DIR}"
for format in png pdf svg; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        pzi_self_highlight.seqids pzi_self_highlight.layout \
        --figsize=18x20 --format="${format}" --notex --keep-chrlabels \
        --chrstyle=roundrect --shadestyle=curve \
        --outfile="pzi_self_highlight_ppa_psp_ppe_pcam_gtr.${format}"
done

echo "[OK] Native JCVI SVG: ${OUT_DIR}/pzi_self_highlight_ppa_psp_ppe_pcam_gtr.svg"
echo "[OK] Native JCVI PDF: ${OUT_DIR}/pzi_self_highlight_ppa_psp_ppe_pcam_gtr.pdf"
echo "[OK] Highlight audit: ${OUT_DIR}/highlight_blocks.tsv"
