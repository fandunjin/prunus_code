#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Combine the uniformly re-estimated CodeML Ks values and draw all five Figure
# 2 panels. This job intentionally does not read the legacy MAFFT/YN outputs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
INPUT_ROOT="${INPUT_ROOT:-${BASE}/04_codeml_inputs}"
RESULT_ROOT="${RESULT_ROOT:-${BASE}/05_codeml_results}"
PLOT_DATA="${PLOT_DATA:-${BASE}/06_codeml_plot_data}"
OUT="${OUT:-${BASE}/07_codeml_five_ks_figure}"
PLOT_PYTHON="${PLOT_PYTHON:-/opt/plant-wakr-venv/bin/python3}"
XMAX="${XMAX:-5.0}"
SELF_PLOT_XMAX="${SELF_PLOT_XMAX:-3.0}"
SELF_AXIS_XMAX="${SELF_AXIS_XMAX:-5.0}"
SELF_MIN_KS="${SELF_MIN_KS:-0.02}"
BANDWIDTH="${BANDWIDTH:-0.035}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${INPUT_ROOT}/COMPLETE"
need "${RESULT_ROOT}/COMPLETE"
need "${ROOT_DIR}/scripts/combine_fig2_repartition_codeml_results.py"
need "${ROOT_DIR}/scripts/plot_figure2_repartition_five_ks.py"
[[ -x "${PLOT_PYTHON}" ]] || die "Python unavailable: ${PLOT_PYTHON}"
"${PLOT_PYTHON}" -c 'import matplotlib, numpy; print("matplotlib=" + matplotlib.__version__)' \
  || die "PLOT_PYTHON must provide matplotlib and numpy"

archive_if_needed() {
  local path="$1"
  if [[ -e "${path}" ]] && [[ -n "$(find "${path}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    if [[ "${FORCE}" != "1" ]]; then
      die "Existing output: ${path}. Set FORCE=1 to archive it before regeneration."
    fi
    mv "${path}" "${path}.replaced.$(date +%Y%m%d_%H%M%S)"
  fi
}

archive_if_needed "${PLOT_DATA}"
archive_if_needed "${OUT}"

"${PLOT_PYTHON}" "${ROOT_DIR}/scripts/combine_fig2_repartition_codeml_results.py" \
  --input-root "${INPUT_ROOT}" --results-root "${RESULT_ROOT}" --outdir "${PLOT_DATA}"
need "${PLOT_DATA}/COMPLETE"
need "${PLOT_DATA}/five_panel_ks.tsv"

"${PLOT_PYTHON}" "${ROOT_DIR}/scripts/plot_figure2_repartition_five_ks.py" \
  --base "${BASE}" --codeml-source "${PLOT_DATA}" --outdir "${OUT}" \
  --xmax "${XMAX}" --self-plot-xmax "${SELF_PLOT_XMAX}" \
  --self-axis-xmax "${SELF_AXIS_XMAX}" --self-min-ks "${SELF_MIN_KS}" \
  --bandwidth "${BANDWIDTH}"

for extension in pdf svg png tiff; do
  need "${OUT}/figures/figure2_repartition_five_ks.${extension}"
done
for panel in a b c d e; do
  for extension in pdf svg png tiff; do
    need "${OUT}/figures/figure2_repartition_ks_${panel}.${extension}"
  done
done
need "${OUT}/source_data/ks_plot_source.tsv"
need "${OUT}/source_data/comparison_audit.tsv"
need "${OUT}/source_data/curve_summary.tsv"
need "${OUT}/source_data/kde_peak_summary.tsv"
grep -q 'PAML CodeML pairwise F3X4' "${OUT}/method_contract.txt" \
  || die "Formal figure method contract is not CodeML"

echo "[OK] Uniform CodeML Figure 2 Ks panels: ${OUT}/figures"
echo "[OK] Pair-level source data and audit: ${OUT}/source_data"
