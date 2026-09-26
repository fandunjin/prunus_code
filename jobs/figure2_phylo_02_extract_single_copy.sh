#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
INPUTS="${BASE}/01_frozen_inputs"
ORTHO="${BASE}/02_orthofinder"
OUT="${BASE}/03_single_copy"
PYTHON_BIN="${PYTHON_BIN:-python3}"

test -e "${ORTHO}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_01 first." >&2; exit 3; }
test ! -e "${OUT}" || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 3; }
mkdir -p "${OUT}"

for dataset in main sensitivity; do
    result_file="${ORTHO}/${dataset}.results_path.txt"
    test -s "${result_file}" || {
        echo "[ERROR] Missing or empty results path: ${result_file}" >&2
        exit 3
    }
    IFS= read -r result < "${result_file}"
    test -d "${result}/Orthogroups" || {
        echo "[ERROR] Invalid OrthoFinder result directory: ${result}" >&2
        exit 3
    }
    "${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_extract_orthogroups.py" \
        --orthogroups "${result}/Orthogroups/Orthogroups.tsv" \
        --pep-dir "${INPUTS}/${dataset}/pep" \
        --out-dir "${OUT}/${dataset}" \
        --min-occupancy 1.0 \
        --min-genes 200 \
        --required pzi_A,pzi_B,psp_A,psp_B,ppa_A,ppa_B,gtr
done
touch "${OUT}/COMPLETE"
