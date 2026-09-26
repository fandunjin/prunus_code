#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Formal Figure 4 A/B fractionation analysis and Python visualization.
# Run inside the local WSL distribution:
#   bash jobs/figure4_ab_fractionation_formal_local.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-/opt/plant-wakr-venv/bin/python3}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure4_ab_fractionation_reference_style}"
SCRIPT="${ROOT_DIR}/scripts/plot_figure4_ab_fractionation_formal.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${SCRIPT}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import numpy, matplotlib' >/dev/null 2>&1 || die "Python needs numpy and matplotlib: ${PYTHON_BIN}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set OUT_ROOT for a rerun"

mkdir -p "${OUT_ROOT}"
"${PYTHON_BIN}" "${SCRIPT}" --root "${ROOT_DIR}" --outdir "${OUT_ROOT}"

printf '%s\n' \
  "backend=Python/matplotlib in local WSL" \
  "analysis=Figure 4a-4e formal A/B fractionation workflow" \
  "output=${OUT_ROOT}" \
  > "${OUT_ROOT}/job_manifest.txt"

echo "[OK] Formal Figure 4 workflow completed: ${OUT_ROOT}"
