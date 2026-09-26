#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Combine formal upper/lower Ks audits and render Figure 2a in Python.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"

die() { echo "[ERROR] $*" >&2; exit 2; }

PYTHON_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
UPPER_ROOT="${BASE}/05_wgd_ks_upper"
LOWER_ROOT="${BASE}/07_wgd_ks_lower"
DATA_OUT="${BASE}/08_panel_a_reference_data"
PLOT_OUT="${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2a_reference_exact"
COMBINER="${ROOT_DIR}/scripts/fig2_panel_a_combine_reference_ks.py"
PLOTTER="${ROOT_DIR}/scripts/plot_fig2_panel_a_reference_ks.py"
FORCE_REPLOT="${FORCE_REPLOT:-0}"

for file in "${COMBINER}" "${PLOTTER}"; do [[ -s "${file}" ]] || die "Missing script: ${file}"; done
for species in PZI PSP PPA; do [[ -e "${UPPER_ROOT}/${species}/COMPLETE" ]] || die "Missing upper completion: ${species}"; done
for track in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B PCAM PPE; do [[ -e "${LOWER_ROOT}/${track}/COMPLETE" ]] || die "Missing lower completion: ${track}"; done
if [[ -e "${PLOT_OUT}/COMPLETE" && "${FORCE_REPLOT}" != "1" ]]; then
  die "Figure output already complete: ${PLOT_OUT}; use FORCE_REPLOT=1 to regenerate"
fi

if [[ -e "${DATA_OUT}/COMPLETE" ]]; then
  echo "[OK] Reusing combined Figure 2a source tables: ${DATA_OUT}"
else
  mkdir -p "${DATA_OUT}"
  "${PYTHON_BIN}" "${COMBINER}" \
    --upper-root "${UPPER_ROOT}" \
    --lower-root "${LOWER_ROOT}" \
    --outdir "${DATA_OUT}"
fi

mkdir -p "${PLOT_OUT}"

PLOT_ARGS=(
  --upper "${DATA_OUT}/panel_a_upper_ks.tsv"
  --lower "${DATA_OUT}/panel_a_lower_ks.tsv"
  --outdir "${PLOT_OUT}"
  --bandwidth 0.04
  --upper-xmax 0.6
  --lower-xmax 1.0
)
if [[ "${FORCE_REPLOT}" == "1" ]]; then
  PLOT_ARGS+=(--overwrite)
fi
MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen "${PYTHON_BIN}" "${PLOTTER}" "${PLOT_ARGS[@]}"

for file in "${PLOT_OUT}/figure2_panel_a_reference_ks.svg" "${PLOT_OUT}/figure2_panel_a_reference_ks.pdf" "${PLOT_OUT}/figure2_panel_a_reference_ks.tiff"; do
  [[ -s "${file}" ]] || die "Missing plot output: ${file}"
done
printf 'Formal Figure 2a combined data and plots completed.\n' > "${PLOT_OUT}/COMPLETE"
echo "[OK] Formal Figure 2a data: ${DATA_OUT}"
echo "[OK] Formal Figure 2a plots: ${PLOT_OUT}"
