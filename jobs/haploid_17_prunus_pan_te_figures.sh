#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
POSTPROCESS_DIR="${POSTPROCESS_DIR:-${OUT_ROOT}/postprocess}"
FIGURE_DIR="${FIGURE_DIR:-${OUT_ROOT}/figures}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/plot_prunus_pan_te.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
MAX_DRIVERS="${MAX_DRIVERS:-18}"
AGE_FAMILIES="${AGE_FAMILIES:-4}"
MAX_DIVERGENCE="${MAX_DIVERGENCE:-50}"
MIN_LTR_OVERLAP_FRACTION="${MIN_LTR_OVERLAP_FRACTION:-0.50}"
MIN_INTACT_ELEMENTS="${MIN_INTACT_ELEMENTS:-20}"
MIN_LTR_LIBRARY_LENGTH="${MIN_LTR_LIBRARY_LENGTH:-500}"
SUBSTITUTION_RATE="${SUBSTITUTION_RATE:-}"

[ -s "${CONDA_SH}" ] || { echo "ERROR: missing conda initialization: ${CONDA_SH}" >&2; exit 1; }
[ -s "${SCRIPT}" ] || { echo "ERROR: missing plotting script: ${SCRIPT}" >&2; exit 1; }
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

args=(
    --postprocess-dir "${POSTPROCESS_DIR}"
    --output-dir "${FIGURE_DIR}"
    --max-drivers "${MAX_DRIVERS}"
    --age-families "${AGE_FAMILIES}"
    --max-divergence "${MAX_DIVERGENCE}"
    --min-ltr-overlap-fraction "${MIN_LTR_OVERLAP_FRACTION}"
    --min-intact-elements "${MIN_INTACT_ELEMENTS}"
    --min-ltr-library-length "${MIN_LTR_LIBRARY_LENGTH}"
)
if [ -n "${SUBSTITUTION_RATE}" ]; then
    args+=(--substitution-rate "${SUBSTITUTION_RATE}")
fi

python "${SCRIPT}" "${args[@]}"

echo "=== Step 17 completed ==="
