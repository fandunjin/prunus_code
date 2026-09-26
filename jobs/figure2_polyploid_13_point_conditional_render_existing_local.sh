#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Render existing completed POInT results without rerunning likelihood fits.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
POINT_ROOT="${POINT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_point_conditional_fixed_r1}"
PARTITION_ROOT="${PARTITION_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/07_point_partition_models_r1}"
OFFICIAL_OUT="${OFFICIAL_OUT:-${POINT_ROOT}/05_panel_d}"
INTEGRATED_OUT="${INTEGRATED_OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/10_point_panel_d_integrated_conditional}"
PLOT_PYTHON_BIN="${PLOT_PYTHON_BIN:-}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${POINT_ROOT}/04_summary/point_model_fit_all.tsv"
require_file "${POINT_ROOT}/04_summary/point_branch_transition_summary.tsv"
require_file "${POINT_ROOT}/04_summary/point_posterior_qc.tsv"
require_file "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv"
require_file "${PARTITION_ROOT}/03_summary/five_model_comparison.tsv"

if [[ -z "${PLOT_PYTHON_BIN}" ]]; then
  for candidate in \
      "/opt/plant-wakr-venv/bin/python3" \
      "${CONDA_ROOT}/envs/wgd2-reference-clean/bin/python" \
      "python3"; do
    if { [[ "${candidate}" == */* && -x "${candidate}" ]] || command -v "${candidate}" >/dev/null 2>&1; } \
        && "${candidate}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
      PLOT_PYTHON_BIN="${candidate}"
      break
    fi
  done
fi
[[ -n "${PLOT_PYTHON_BIN}" ]] || die "No Python interpreter with matplotlib and numpy found; set PLOT_PYTHON_BIN"

[[ ! -e "${OFFICIAL_OUT}" ]] || die "Output directory already exists: ${OFFICIAL_OUT}; set OFFICIAL_OUT for a rerun"
[[ ! -e "${INTEGRATED_OUT}" ]] || die "Output directory already exists: ${INTEGRATED_OUT}; set INTEGRATED_OUT for a rerun"

mkdir -p "${OFFICIAL_OUT}"
"${PLOT_PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_panel_d_point_official.py" \
  --fit "${POINT_ROOT}/04_summary/point_model_fit_all.tsv" \
  --transition "${POINT_ROOT}/04_summary/point_branch_transition_summary.tsv" \
  --posterior-qc "${POINT_ROOT}/04_summary/point_posterior_qc.tsv" \
  --posterior-summary "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv" \
  --outdir "${OFFICIAL_OUT}"

POINT_ROOT="${POINT_ROOT}" \
PARTITION_ROOT="${PARTITION_ROOT}" \
OUT_ROOT="${INTEGRATED_OUT}" \
PLOT_PYTHON_BIN="${PLOT_PYTHON_BIN}" \
  bash "${ROOT_DIR}/jobs/figure2_polyploid_11_point_panel_d_integrated_local.sh"

echo "[OK] Rendered completed POInT conditional results: ${INTEGRATED_OUT}"
