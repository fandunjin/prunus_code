#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
INTEGRATE="${ROOT_DIR}/scripts/prunus_pan_te_tandem_integrated_summary.py"
PLOT="${ROOT_DIR}/scripts/plot_prunus_pan_te_tandem_te_integrated.py"
TE_DIR="${PAN_ROOT}/postprocess_v2"
TANDEM_DIR="${PAN_ROOT}/figure1_tandem"
OUT="${PAN_ROOT}/figure1_tandem_te_integrated"
FIG_OUT="${PAN_ROOT}/figures_tandem_te_integrated"

for path in "${PYTHON}" "${CONDA_SH}" "${INTEGRATE}" "${PLOT}" "${TE_DIR}/te_composition.tsv" "${TANDEM_DIR}/main/tandem_main_partition.tsv" "${TANDEM_DIR}/score50/tandem_period_summary.tsv"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
mkdir -p "${OUT}" "${FIG_OUT}"
"${PYTHON}" "${INTEGRATE}" \
    --te-composition "${TE_DIR}/te_composition.tsv" \
    --tandem-partition "${TANDEM_DIR}/main/tandem_main_partition.tsv" \
    --tandem-period-summary "${TANDEM_DIR}/score50/tandem_period_summary.tsv" \
    --output-dir "${OUT}"
source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
export MPLBACKEND=Agg
python -c 'import matplotlib, numpy; print("matplotlib", matplotlib.__version__); print("numpy", numpy.__version__)'
python "${PLOT}" \
    --integrated-dir "${OUT}" \
    --te-composition "${TE_DIR}/te_composition.tsv" \
    --output-dir "${FIG_OUT}"
echo "=== Integrated pan-TE/tandem Figure 1 export completed ==="
