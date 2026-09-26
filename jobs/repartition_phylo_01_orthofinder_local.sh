#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
INPUTS="${BASE}/01_frozen_inputs"
OUT_ROOT="${OUT_ROOT:-${BASE}/02_orthofinder}"
ORTHOFINDER_BIN="${ORTHOFINDER_BIN:-orthofinder}"
THREADS="${THREADS:-12}"
ANALYSIS_THREADS="${ANALYSIS_THREADS:-4}"

[[ -e "${INPUTS}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_00_freeze_inputs_local.sh first" >&2; exit 2; }
[[ ! -e "${OUT_ROOT}" ]] || { echo "[ERROR] Output exists: ${OUT_ROOT}" >&2; exit 2; }
command -v "${ORTHOFINDER_BIN}" >/dev/null 2>&1 || { echo "[ERROR] OrthoFinder unavailable: ${ORTHOFINDER_BIN}" >&2; exit 2; }
mkdir -p "${OUT_ROOT}"

for dataset in main sensitivity; do
  "${ORTHOFINDER_BIN}" -f "${INPUTS}/${dataset}/pep" -S diamond -t "${THREADS}" -a "${ANALYSIS_THREADS}" -o "${OUT_ROOT}/${dataset}"
  result=""
  for candidate in "${OUT_ROOT}/${dataset}" "${OUT_ROOT}/${dataset}"/Results* "${OUT_ROOT}/${dataset}"/OrthoFinder/Results*; do
    if [[ -s "${candidate}/Orthogroups/Orthogroups.tsv" || -s "${candidate}/Orthogroups/Orthogroups.txt" ]]; then result="${candidate}"; break; fi
  done
  [[ -n "${result}" ]] || { echo "[ERROR] OrthoFinder result not found for ${dataset}" >&2; exit 2; }
  printf '%s\n' "${result}" > "${OUT_ROOT}/${dataset}.results_path.txt"
done
{
  echo "orthofinder=${ORTHOFINDER_BIN}"
  "${ORTHOFINDER_BIN}" -h 2>&1 | head -n 3 || true
  echo "threads=${THREADS}"
  echo "analysis_threads=${ANALYSIS_THREADS}"
  echo "main_input=${INPUTS}/main/pep"
  echo "sensitivity_input=${INPUTS}/sensitivity/pep"
} > "${OUT_ROOT}/software_versions.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] New repartitioned OrthoFinder runs completed: ${OUT_ROOT}"
