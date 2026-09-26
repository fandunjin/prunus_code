#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
FIGURE_DIR="${FIGURE_DIR:-${OUT_ROOT}/figures_followup}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/plot_prunus_pan_te_followup.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"

for path in "${SCRIPT}" "${CONDA_SH}" \
    "${FOLLOWUP_ROOT}/component_audit/component_audit_summary.tsv" \
    "${FOLLOWUP_ROOT}/sololtr_summary/sololtr_by_target_lineage.tsv" \
    "${FOLLOWUP_ROOT}/orthosite/combined/orthosite_family_summary.tsv"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
export MPLBACKEND=Agg

python - <<'PY'
import matplotlib
import numpy
import PIL
print("matplotlib", matplotlib.__version__)
print("numpy", numpy.__version__)
print("Pillow", PIL.__version__)
PY

python "${SCRIPT}" \
    --followup-root "${FOLLOWUP_ROOT}" \
    --output-dir "${FIGURE_DIR}"

echo "=== Follow-up TE figure export completed ==="
