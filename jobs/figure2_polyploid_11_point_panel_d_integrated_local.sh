#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Render the final integrated POInT evidence panel for Figure 2d.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
POINT_ROOT="${POINT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_official_reference_exact_r13}"
PARTITION_ROOT="${PARTITION_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/07_point_partition_models_r1}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/08_point_panel_d_integrated}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PLOT_PYTHON_BIN="${PLOT_PYTHON_BIN:-}"
PLOT="${ROOT_DIR}/scripts/plot_fig2_panel_d_point_integrated.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${POINT_ROOT}/04_summary/point_model_fit_all.tsv"
require_file "${PARTITION_ROOT}/03_summary/five_model_comparison.tsv"
require_file "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv"
[[ -f "${PLOT}" ]] || die "Missing plot script: ${PLOT}"

if [[ -z "${PLOT_PYTHON_BIN}" ]]; then
  for candidate in "${PYTHON_BIN}" \
      "/opt/plant-wakr-venv/bin/python3" \
      "${CONDA_ROOT}/envs/wgd2-reference-clean/bin/python"; do
    if [[ -x "${candidate}" ]] && "${candidate}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
      PLOT_PYTHON_BIN="${candidate}"
      break
    fi
  done
fi
[[ -n "${PLOT_PYTHON_BIN}" ]] || die "No Python interpreter with matplotlib and numpy found; set PLOT_PYTHON_BIN"
if [[ -e "${OUT_ROOT}" ]] && [[ -n "$(find "${OUT_ROOT}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  die "Output directory is not empty: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"
fi

mkdir -p "${OUT_ROOT}"
TRANSITION_FILE="${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_branch_transition_summary.tsv"
POSTERIOR_PILLARS=$(awk 'NR > 1 && NF {n++} END {print n + 0}' "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv")
transition_args=()
conditional_status="unavailable"
if [[ -s "${TRANSITION_FILE}" ]] && awk 'NR > 1 && NF { found=1 } END { exit !found }' "${TRANSITION_FILE}"; then
  transition_args=(--transition "${TRANSITION_FILE}")
  conditional_status="available"
fi
"${PLOT_PYTHON_BIN}" "${PLOT}" \
  --fit "${POINT_ROOT}/04_summary/point_model_fit_all.tsv" \
  --partitions "${PARTITION_ROOT}/03_summary/five_model_comparison.tsv" \
  --posterior "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv" \
  "${transition_args[@]}" \
  --outdir "${OUT_ROOT}"

cat > "${OUT_ROOT}/software_and_parameters.txt" <<EOF
method=official POInT v1.63
fit=${POINT_ROOT}/04_summary/point_model_fit_all.tsv
partition_models=${PARTITION_ROOT}/03_summary/five_model_comparison.tsv
posterior=${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv
posterior_pillars=${POSTERIOR_PILLARS}
panel_1=Delta AIC, symmetric versus WGD_bias_fix
panel_2=Delta AIC bars and Delta BIC points, five event partitions
panel_3=conditional branch-transition expectations when available
panel_4=best posterior tracking-state fractions
conditional_output=${conditional_status}
conditional_transition_summary=${TRANSITION_FILE}
EOF
echo "[OK] Integrated Figure 2d completed: ${OUT_ROOT}"
