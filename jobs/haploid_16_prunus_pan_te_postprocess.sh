#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
SUMMARY_DIR="${SUMMARY_DIR:-${OUT_ROOT}/summary}"
DISCOVERY_DIR="${DISCOVERY_DIR:-${OUT_ROOT}/discovery}"
PAN_LIBRARY="${PAN_LIBRARY:-${OUT_ROOT}/library/frozen/prunus_pan_te_v1.fa}"
POSTPROCESS_DIR="${POSTPROCESS_DIR:-${OUT_ROOT}/postprocess}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_postprocess.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
THREADS="${THREADS:-4}"
CUMULATIVE_THRESHOLD="${CUMULATIVE_THRESHOLD:-0.80}"
DIVERGENCE_BIN_WIDTH="${DIVERGENCE_BIN_WIDTH:-1.0}"
# Leave empty until a defensible Prunus neutral substitution rate is fixed.
SUBSTITUTION_RATE="${SUBSTITUTION_RATE:-}"

[ -s "${CONDA_SH}" ] || { echo "ERROR: missing conda initialization: ${CONDA_SH}" >&2; exit 1; }
[ -s "${SCRIPT}" ] || { echo "ERROR: missing post-process script: ${SCRIPT}" >&2; exit 1; }
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"

args=(
    --summary-dir "${SUMMARY_DIR}"
    --discovery-dir "${DISCOVERY_DIR}"
    --pan-library "${PAN_LIBRARY}"
    --output-dir "${POSTPROCESS_DIR}"
    --threads "${THREADS}"
    --cumulative-threshold "${CUMULATIVE_THRESHOLD}"
    --divergence-bin-width "${DIVERGENCE_BIN_WIDTH}"
)
if [ -n "${SUBSTITUTION_RATE}" ]; then
    args+=(--substitution-rate "${SUBSTITUTION_RATE}")
fi

python "${SCRIPT}" "${args[@]}"

echo "=== Step 16 completed ==="
