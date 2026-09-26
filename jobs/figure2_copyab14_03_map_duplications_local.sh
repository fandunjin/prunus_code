#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RESULT_DIR="${RESULT_DIR:-}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/0903/analysisphylogenetics/03_copyAB_duplication_mapping}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ -z "${RESULT_DIR}" && -s "${ROOT_DIR}/0903/analysisphylogenetics/02_orthofinder_copyAB14/latest.txt" ]]; then RESULT_DIR="$(head -n 1 "${ROOT_DIR}/0903/analysisphylogenetics/02_orthofinder_copyAB14/latest.txt")"; fi
[[ -n "${RESULT_DIR}" && -d "${RESULT_DIR}" ]] || { echo "[ERROR] Set RESULT_DIR to OrthoFinder result directory" >&2; exit 2; }
mkdir -p "${OUT_DIR}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/summarize_copyab_duplications.py" --result "${RESULT_DIR}" --outdir "${OUT_DIR}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_copyab_duplication_mapping.py" --tree "${RESULT_DIR}/Species_Tree/SpeciesTree_rooted_node_labels.txt" --counts "${OUT_DIR}/copyAB_duplication_node_counts.tsv" --outdir "${OUT_DIR}"
