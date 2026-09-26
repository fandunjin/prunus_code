#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
SCRIPT="${ROOT_DIR}/scripts/plot_prunus_pan_te_tandem_main.py"
OUT="${PAN_ROOT}/figures_tandem_main"

for path in "${CONDA_SH}" "${SCRIPT}" "${FIG_ROOT}/main/tandem_main_partition.tsv" "${FIG_ROOT}/score50/tandem_period_summary.tsv" "${FIG_ROOT}/loci/pzi_tandem_locus_summary.tsv" "${FIG_ROOT}/main/tandem_score_sensitivity.tsv"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
for target in pziA1 pziA2 pziB1 pziB2; do
    [ -s "${FIG_ROOT}/ont/${target}.ont_support.tsv" ] || { echo "ERROR: missing ONT support table for ${target}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
export MPLBACKEND=Agg
python -c 'import matplotlib, numpy; print("matplotlib", matplotlib.__version__); print("numpy", numpy.__version__)'
python "${SCRIPT}" \
    --main-dir "${FIG_ROOT}/main" \
    --tandem-summary-dir "${FIG_ROOT}/score50" \
    --locus-dir "${FIG_ROOT}/loci" \
    --ont-dir "${FIG_ROOT}/ont" \
    --sensitivity "${FIG_ROOT}/main/tandem_score_sensitivity.tsv" \
    --output-dir "${OUT}"
echo "=== Tandem main-figure exports completed ==="
