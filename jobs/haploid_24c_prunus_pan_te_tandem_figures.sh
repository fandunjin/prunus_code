#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/plot_prunus_pan_te_tandem.py}"
SUMMARY_DIR="${SUMMARY_DIR:-${PAN_ROOT}/followup/tandem_summary}"
OUT="${OUT:-${PAN_ROOT}/figures_tandem}"
[ -s "${CONDA_SH}" ] && [ -s "${SCRIPT}" ] || { echo "ERROR: missing conda/script" >&2; exit 1; }
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
export MPLBACKEND=Agg
python -c 'import matplotlib, numpy; print("matplotlib", matplotlib.__version__); print("numpy", numpy.__version__)'
python "${SCRIPT}" --summary-dir "${SUMMARY_DIR}" --output-dir "${OUT}" --top-motifs 12
echo "=== Tandem-repeat figure export completed ==="
