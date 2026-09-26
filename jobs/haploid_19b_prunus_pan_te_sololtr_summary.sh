#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_sololtr_summary.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"

for path in "${SCRIPT}" "${CONDA_SH}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"

python "${SCRIPT}" \
    --sololtr-root "${FOLLOWUP_ROOT}/sololtr" \
    --output-dir "${FOLLOWUP_ROOT}/sololtr_summary"

echo "=== soloLTR summary completed ==="
