#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Replot all nine GTR chromosomes while retaining every PZI-A/PZI-B chromosome.
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-/opt/plant-wakr-venv/bin/python3}"
FORMAL_ROOT="${FORMAL_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact}"
OUT_ROOT="${OUT_ROOT:-${FORMAL_ROOT}/07_all_gtr_chromosome_audit_all_pzi}"
WINDOW_ROOT="${OUT_ROOT}/01_windows"
PLOT_ROOT="${OUT_ROOT}/02_plot"
ANCHORS="${ANCHORS:-${FORMAL_ROOT}/04_quota_align/normalized/quota_filtered_anchors.tsv}"
CATALOG="${CATALOG:-${FORMAL_ROOT}/01_catalog}"

die(){ echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${ANCHORS}" ]] || die "Missing formal quota-filtered anchors: ${ANCHORS}"
[[ -s "${CATALOG}/gtr_gene_order.tsv" ]] || die "Missing GTR catalog"
[[ -s "${CATALOG}/query_gene_track_map.tsv" ]] || die "Missing query catalog"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}"
mkdir -p "${OUT_ROOT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_fractbias_all_gtr_all_pzi.py" \
  --anchors "${ANCHORS}" --catalog "${CATALOG}" --outdir "${WINDOW_ROOT}" --window-size 100
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_panel_c_all_gtr_all_pzi.py" \
  --windows "${WINDOW_ROOT}/all_gtr_all_pzi_fractbias_windows.tsv" \
  --selection "${WINDOW_ROOT}/all_gtr_all_pzi_selected_tracks.tsv" \
  --outdir "${PLOT_ROOT}"
printf '[OK] All-PZI 9-panel and individual GTR retention plots completed: %s\n' "${PLOT_ROOT}"
