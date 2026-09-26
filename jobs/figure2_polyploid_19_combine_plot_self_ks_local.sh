#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Combine six genome-native self-syntenic paralog Ks analyses, then render a
# compact two-row diagnostic figure.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
FORCE_REPLOT="${FORCE_REPLOT:-0}"
PYTHON_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin/python3"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"

SELF_BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_within_genome_ks"
DATA_OUT="${SELF_BASE}/06_combined_source_data"
PLOT_OUT="${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2_within_genome_ks"
COMBINER="${ROOT_DIR}/scripts/fig2_combine_within_genome_ks.py"
PLOTTER="${ROOT_DIR}/scripts/plot_fig2_within_genome_ks.py"

for file in "${COMBINER}" "${PLOTTER}"; do
  [[ -s "${file}" ]] || die "Missing input: ${file}"
done
for label in PZI PSP PPA PPE PCAM GTR; do
  [[ -e "${SELF_BASE}/05_wgd_ks/${label}/COMPLETE" ]] || die "Missing completed self-syntenic Ks: ${label}"
done

if [[ ! -e "${DATA_OUT}/COMPLETE" ]]; then
  mkdir -p "${DATA_OUT}"
  "${PYTHON_BIN}" "${COMBINER}" \
    --self-root "${SELF_BASE}/05_wgd_ks" \
    --outdir "${DATA_OUT}"
else
  echo "[INFO] Reusing combined source data: ${DATA_OUT}"
fi
if [[ -e "${PLOT_OUT}/COMPLETE" && "${FORCE_REPLOT}" != "1" ]]; then
  die "Plot output already complete: ${PLOT_OUT}; use FORCE_REPLOT=1 to regenerate"
fi
mkdir -p "${PLOT_OUT}"
PLOT_ARGS=(--source "${DATA_OUT}/within_genome_ks.tsv" --outdir "${PLOT_OUT}" --bandwidth 0.035 --xmax 3.0)
if [[ "${FORCE_REPLOT}" == "1" ]]; then PLOT_ARGS+=(--overwrite); fi
MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen "${PYTHON_BIN}" "${PLOTTER}" "${PLOT_ARGS[@]}"
for file in "${PLOT_OUT}/figure2_within_genome_ks.svg" "${PLOT_OUT}/figure2_within_genome_ks.pdf" "${PLOT_OUT}/figure2_within_genome_ks.tiff"; do
  [[ -s "${file}" ]] || die "Missing plot output: ${file}"
done
printf 'Within-genome Ks diagnostic figure completed.\n' > "${PLOT_OUT}/COMPLETE"
echo "[OK] Source data: ${DATA_OUT}"
echo "[OK] Figure: ${PLOT_OUT}"
