#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Plot native OrthoFinder duplication mappings on its A/B-labelled species tree.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OF_ROOT="${OF_ROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13}"
OUT="${OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2b_orthofinder_native_duplications}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

TREE="${OF_ROOT}/Species_Tree/SpeciesTree_rooted_node_labels.txt"
COUNTS="${OF_ROOT}/Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv"
[[ -s "${TREE}" ]] || { echo "[ERROR] Missing OrthoFinder labelled tree: ${TREE}" >&2; exit 2; }
[[ -s "${COUNTS}" ]] || { echo "[ERROR] Missing OrthoFinder duplication summary: ${COUNTS}" >&2; exit 2; }

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_orthofinder_duplication_mapping.py" \
    --tree "${TREE}" \
    --counts "${COUNTS}" \
    --outdir "${OUT}"
