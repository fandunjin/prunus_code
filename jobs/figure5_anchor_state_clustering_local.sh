#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-/opt/plant-wakr-venv/bin/python3}"
PILLARS="${PILLARS:-${ROOT_DIR}/downstream_analysis/figure4_ab_fractionation_reference_style_v4/panel4a_anchor_state_by_pillar.tsv}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure5_anchor_state_clustering}"
PERMUTATIONS="${PERMUTATIONS:-2000}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${PILLARS}" ]] || die "Missing strict POInT pillar table: ${PILLARS}"
[[ -s "${ROOT_DIR}/scripts/analyze_anchor_state_clustering.py" ]] || die "Missing analysis script"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}"

mkdir -p "${OUT_ROOT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/analyze_anchor_state_clustering.py" \
  --pillars "${PILLARS}" \
  --outdir "${OUT_ROOT}" \
  --permutations "${PERMUTATIONS}" \
  --seed 20260828

printf '[OK] Figure 5 anchor-state clustering completed: %s\n' "${OUT_ROOT}"
