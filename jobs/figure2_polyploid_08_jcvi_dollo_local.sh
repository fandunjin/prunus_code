#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for the JCVI-anchored Dollo-parsimony Fig. 2d analysis.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STATE_TABLE="${STATE_TABLE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/alltaxa_fractionation_states/summary/gtr_fractionation_states.tsv}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_jcvi_dollo_fractionation_v2}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CALC_SCRIPT="${ROOT_DIR}/scripts/fig2_jcvi_dollo_fractionation.py"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_fig2_panel_d_point_reference_style.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input/output: $1"; }

for file in "${STATE_TABLE}" "${CALC_SCRIPT}" "${PLOT_SCRIPT}"; do require_file "${file}"; done
"${PYTHON_BIN}" -c 'import matplotlib' >/dev/null 2>&1 || die "matplotlib unavailable in ${PYTHON_BIN}"
mkdir -p "${OUT_ROOT}"
if [[ ! -e "${OUT_ROOT}/01_dollo/COMPLETE" ]]; then
    "${PYTHON_BIN}" "${CALC_SCRIPT}" --states "${STATE_TABLE}" --outdir "${OUT_ROOT}/01_dollo"
fi
[[ ! -e "${OUT_ROOT}/02_panel_d" ]] || die "Plot output already exists: ${OUT_ROOT}/02_panel_d"
"${PYTHON_BIN}" "${PLOT_SCRIPT}" \
    --point-output "${OUT_ROOT}/01_dollo/jcvi_dollo_branch_summary.tsv" \
    --mode jcvi_dollo --outdir "${OUT_ROOT}/02_panel_d"
{
    echo "method=JCVI conservative pillars plus Dollo parsimony"
    echo "tree=(GTR,(PPA,(PSP,PZI)))"
    echo "source=${STATE_TABLE}"
    echo "boundary=anchor attrition, not POInT likelihood-based gene loss"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] JCVI-Dollo Fig. 2d completed: ${OUT_ROOT}"
