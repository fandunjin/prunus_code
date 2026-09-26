#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Python-only formal panel b plot from i-ADHoRe-derived PUG results.

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
PUG_ROOT="${PUG_ROOT:-${BASE}/03_pug_reference_iadhore_allastral}"
OUT="${OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2b_pug_reference_iadhore_allastral}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[[ -e "${PUG_ROOT}/COMPLETE" ]] || { echo "[ERROR] Missing PUG completion: ${PUG_ROOT}/COMPLETE" >&2; exit 2; }
[[ -s "${PUG_ROOT}/pug_node_counts.tsv" ]] || { echo "[ERROR] Missing PUG summary" >&2; exit 2; }
MAPPED="${PUG_ROOT}/shared_candidate/fig2b_shared_candidate_Paralog_Pairs_Nodes_Bootstraps.txt"
TREE="${PUG_ROOT}/shared_candidate/fig2b_shared_candidate_Labeled_Species_Tree.tre"
[[ -s "${MAPPED}" ]] || { echo "[ERROR] Missing raw pooled PUG rows" >&2; exit 2; }
[[ -s "${TREE}" ]] || { echo "[ERROR] Missing labeled PUG species tree" >&2; exit 2; }
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output already exists: ${OUT}" >&2; exit 2; }

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_panel_b_pug_reference_style_allastral.py" \
    --pug "${PUG_ROOT}/pug_node_counts.tsv" --mapped-rows "${MAPPED}" \
    --species-tree "${TREE}" --outdir "${OUT}"
touch "${OUT}/COMPLETE"
echo "[OK] Formal i-ADHoRe/PUG panel b plotted: ${OUT}"
