#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_orthosite.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"

for path in "${SCRIPT}" "${CONDA_SH}" \
    "${FOLLOWUP_ROOT}/orthosite/pziA1/orthosite_site_evidence.tsv" \
    "${FOLLOWUP_ROOT}/orthosite/pziB1/orthosite_site_evidence.tsv"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"

python "${SCRIPT}" combine \
    --source-dir "${FOLLOWUP_ROOT}/orthosite/pziA1" \
    --source-dir "${FOLLOWUP_ROOT}/orthosite/pziB1" \
    --output-dir "${FOLLOWUP_ROOT}/orthosite/combined"

echo "=== A/B orthosite results combined ==="
