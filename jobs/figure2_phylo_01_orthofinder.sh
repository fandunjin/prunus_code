#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
INPUTS="${BASE}/01_frozen_inputs"
OUT_ROOT="${BASE}/02_orthofinder"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
ORTHO_ENV="${ORTHO_ENV:-orthofinder}"
THREADS="${THREADS:-48}"
ANALYSIS_THREADS="${ANALYSIS_THREADS:-8}"

test -e "${INPUTS}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_00 first." >&2; exit 3; }
test ! -e "${OUT_ROOT}" || { echo "[ERROR] Output exists: ${OUT_ROOT}" >&2; exit 3; }
mkdir -p "${OUT_ROOT}"

set +u
source "${CONDA_SH}"
conda activate "${ORTHO_ENV}"
set -u
command -v orthofinder >/dev/null || { echo "[ERROR] orthofinder unavailable" >&2; exit 2; }
command -v diamond >/dev/null || { echo "[ERROR] diamond unavailable" >&2; exit 2; }

for dataset in main sensitivity; do
    orthofinder \
        -f "${INPUTS}/${dataset}/pep" \
        -S diamond \
        -t "${THREADS}" \
        -a "${ANALYSIS_THREADS}" \
        -o "${OUT_ROOT}/${dataset}"

    result=""
    for candidate in \
        "${OUT_ROOT}/${dataset}" \
        "${OUT_ROOT}/${dataset}"/Results* \
        "${OUT_ROOT}/${dataset}"/OrthoFinder/Results*
    do
        if [ -s "${candidate}/Orthogroups/Orthogroups.tsv" ]; then
            result="${candidate}"
            break
        fi
    done
    test -n "${result}" || { echo "[ERROR] OrthoFinder result not found for ${dataset}" >&2; exit 3; }
    printf '%s\n' "${result}" > "${OUT_ROOT}/${dataset}.results_path.txt"
done

{
    echo "OrthoFinder command: $(command -v orthofinder)"
    orthofinder -h 2>&1 | head -n 3 || true
    echo
    echo "DIAMOND:"
    diamond version 2>&1 || true
} > "${OUT_ROOT}/software_versions.txt"
touch "${OUT_ROOT}/COMPLETE"
