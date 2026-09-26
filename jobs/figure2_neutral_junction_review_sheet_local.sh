#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare, but never automatically approve, the Figure 2 junction review sheet.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BLOCK_ROOT="${BLOCK_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_block_ancestry}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
HIFI_TABLE="${HIFI_TABLE:-}"
OUT="${OUT:-${BLOCK_ROOT}/05_junction_manual_review.tsv}"

[[ -s "${BLOCK_ROOT}/04_candidate_validation/junction_validation_summary.tsv" ]] || {
    echo "[ERROR] Missing stage-4 validation summary." >&2; exit 2;
}
args=(--validation "${BLOCK_ROOT}/04_candidate_validation/junction_validation_summary.tsv" --out "${OUT}")
if [[ -n "${HIFI_TABLE}" ]]; then
    [[ -s "${HIFI_TABLE}" ]] || { echo "[ERROR] Missing HIFI_TABLE: ${HIFI_TABLE}" >&2; exit 2; }
    args+=(--hifi "${HIFI_TABLE}")
fi
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_junction_review_sheet.py" "${args[@]}"
echo "[NEXT] Complete Hi-C coordinate checks and manual review, then set final_status=confirmed_rejoined_junction only for confirmed rows."
