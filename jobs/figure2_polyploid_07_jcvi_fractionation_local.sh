#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for the registered JCVI-adapted Fig. 2c analysis.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PAIR_TABLE="${PAIR_TABLE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/alltaxa_fractionation_states/summary/gtr_alltaxa_syntenic_pairs.tsv}"
GTR_BED="${GTR_BED:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs/gtr/gtr.bed}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_jcvi_fractionation_8track}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CALC_SCRIPT="${ROOT_DIR}/scripts/fig2_jcvi_fractbias_equivalent.py"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_fig2_panel_c_fractbias_reference_style.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input/output: $1"; }

for file in "${PAIR_TABLE}" "${GTR_BED}" "${CALC_SCRIPT}" "${PLOT_SCRIPT}"; do require_file "${file}"; done
"${PYTHON_BIN}" -c 'import matplotlib, numpy' >/dev/null 2>&1 || die "matplotlib/numpy unavailable in ${PYTHON_BIN}"
mkdir -p "${OUT_ROOT}"
if [[ ! -e "${OUT_ROOT}/01_windows/COMPLETE" ]]; then
    "${PYTHON_BIN}" "${CALC_SCRIPT}" --pairs "${PAIR_TABLE}" --gtr-bed "${GTR_BED}" --window 100 --step 5 --outdir "${OUT_ROOT}/01_windows"
fi
[[ ! -e "${OUT_ROOT}/02_panel_c" ]] || die "Plot output already exists: ${OUT_ROOT}/02_panel_c"
"${PYTHON_BIN}" "${PLOT_SCRIPT}" \
    --fractbias "${OUT_ROOT}/01_windows/jcvi_fractbias_windows.tsv" \
    --mode jcvi_adapted --outdir "${OUT_ROOT}/02_panel_c"
{
    echo "method=JCVI screened anchors plus centered 100-gene rolling window"
    echo "source=${PAIR_TABLE}"
    echo "display_pairs=C06_A1/B1,C02_A1/B1,C03_A1/B1"
    echo "tracks=PZI_A,PZI_B,PSP_A,PSP_B,PPA_A,PPA_B,PCAM,PPE"
    echo "boundary=not a FractBias software run; unanchored loci are not verified losses"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] JCVI-adapted Fig. 2c completed: ${OUT_ROOT}"
