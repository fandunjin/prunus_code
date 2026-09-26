#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare and plot Figure 2 Ks panels on continuous smooth nonlinear display axes.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
PLOT_ENV="${PLOT_ENV:-wgd2-reference-clean}"
PYTHON_BIN="${PYTHON_BIN:-${CONDA_BASE}/envs/${PLOT_ENV}/bin/python3}"
FORCE_REPLOT="${FORCE_REPLOT:-0}"
PANEL_C_MIN_KS="${PANEL_C_MIN_KS:-0.02}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

[[ -x "${PYTHON_BIN}" ]] || die "Missing plotting Python: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import matplotlib, numpy; print("[OK] matplotlib=" + matplotlib.__version__)'

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid"
DATA_OUT="${BASE}/strict_wgd_mapping/10_smooth_ks_inputs"
PLOT_OUT="${BASE}/fig2_ks_smooth_axis_compact"
PREP="${ROOT_DIR}/scripts/prepare_fig2_piecewise_ks_inputs.py"
PLOTTER="${ROOT_DIR}/scripts/plot_fig2_piecewise_ks.py"
AB="${BASE}/strict_wgd_mapping/06_panel_a_reference_exact/08_panel_a_reference_data/panel_a_upper_ks.tsv"
REFERENCE="${BASE}/subgenome_reference_ks_YN/02_combined_source_data/subgenome_reference_ks.tsv"
SELF_ROOT="${BASE}/strict_wgd_mapping/09_within_genome_ks/05_wgd_ks"

for file in "${PREP}" "${PLOTTER}" "${AB}" "${REFERENCE}"; do require_file "${file}"; done
for label in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B PCAM PPE GTR; do
  require_file "${SELF_ROOT}/${label}/audit/${label}_SELF.all_pair_ks_audit.tsv"
  if [[ "${label}" =~ ^(PZI|PSP|PPA)_(A|B)$ ]]; then
    require_file "${SELF_ROOT}/${label}/run_parameters.tsv"
    awk -F '\t' '$1 == "analysis_scope" && $2 == "independent_subgenome_self_synteny" { found=1 } END { exit !found }' "${SELF_ROOT}/${label}/run_parameters.tsv" \
      || die "${label} is not an independently run subgenome self-synteny analysis"
  fi
done

if [[ -e "${DATA_OUT}/COMPLETE" ]] \
  && grep -q 'independently run PZI_A' "${DATA_OUT}/method_contract.txt"; then
  echo "[INFO] Reusing prepared smooth-axis Ks inputs: ${DATA_OUT}"
else
  [[ ! -e "${DATA_OUT}" || -z "$(find "${DATA_OUT}" -mindepth 1 -maxdepth 1 -print -quit)" ]] \
    || die "Prepared smooth-axis input directory is incomplete or incompatible: ${DATA_OUT}; preserve it for audit and choose a new DATA_OUT path after inspection"
  "${PYTHON_BIN}" "${PREP}" --ab "${AB}" --reference "${REFERENCE}" --self-root "${SELF_ROOT}" --outdir "${DATA_OUT}"
fi

if [[ -e "${PLOT_OUT}/COMPLETE" && "${FORCE_REPLOT}" != "1" ]]; then
  die "Output is already complete: ${PLOT_OUT}; rerun with FORCE_REPLOT=1"
fi
mkdir -p "${PLOT_OUT}"
PLOT_ARGS=(
  --panel-a "${DATA_OUT}/panel_a_ab_homoeologue_ks.tsv"
  --panel-b "${DATA_OUT}/panel_b_subgenome_reference_ks.tsv"
  --panel-c "${DATA_OUT}/panel_c_subgenome_self_ks.tsv"
  --outdir "${PLOT_OUT}"
  --panel-c-min-ks "${PANEL_C_MIN_KS}"
)
if [[ "${FORCE_REPLOT}" == "1" ]]; then PLOT_ARGS+=(--overwrite); fi
MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen "${PYTHON_BIN}" "${PLOTTER}" "${PLOT_ARGS[@]}"

for extension in pdf svg png tiff; do
  require_file "${PLOT_OUT}/figure2_ks_smooth_axis.${extension}"
done
echo "[OK] Compact Figure 2 smooth-axis Ks PDF: ${PLOT_OUT}/figure2_ks_smooth_axis.pdf"
