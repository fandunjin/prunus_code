#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
INPUTS="${BASE}/01_frozen_inputs"
ORTHO="${BASE}/02_orthofinder"
OUT="${OUT:-${BASE}/03_single_copy}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[[ -e "${ORTHO}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_01_orthofinder_local.sh first" >&2; exit 2; }
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 2; }
mkdir -p "${OUT}"
for dataset in main sensitivity; do
  result_file="${ORTHO}/${dataset}.results_path.txt"
  [[ -s "${result_file}" ]] || { echo "[ERROR] Missing results path: ${result_file}" >&2; exit 2; }
  IFS= read -r result < "${result_file}"
  [[ -s "${result}/Orthogroups/Orthogroups.tsv" || -s "${result}/Orthogroups/Orthogroups.txt" ]] || { echo "[ERROR] Invalid OrthoFinder result: ${result}" >&2; exit 2; }
  if [[ -s "${result}/Orthogroups/Orthogroups.tsv" ]]; then
    orthogroups_file="${result}/Orthogroups/Orthogroups.tsv"
  else
    orthogroups_file="${result}/Orthogroups/Orthogroups.txt"
  fi
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_extract_orthogroups.py" \
    --orthogroups "${orthogroups_file}" \
    --pep-dir "${INPUTS}/${dataset}/pep" \
    --out-dir "${OUT}/${dataset}" \
    --min-occupancy 1.0 --min-genes 200 \
    --required pzi_A,pzi_B,psp_A,psp_B,ppa_A,ppa_B,gtr
done
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned single-copy extraction completed: ${OUT}"
