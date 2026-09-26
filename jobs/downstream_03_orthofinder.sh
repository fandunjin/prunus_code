#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
ORTHO_ENV="${ORTHO_ENV:-orthofinder}"
THREADS="${THREADS:-48}"
ANALYSIS_THREADS="${ANALYSIS_THREADS:-8}"
SEARCH_PROGRAM="${SEARCH_PROGRAM:-diamond}"
RUN_ID="${RUN_ID:-${JSUB_JOBID:-${LSB_JOBID:-$$}}}"

mkdir -p "${ANALYSIS_DIR}/orthofinder/results" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${ORTHO_ENV}"
set -u

for exe in orthofinder "${SEARCH_PROGRAM}"; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in ORTHO_ENV=${ORTHO_ENV}: ${exe}" >&2
        exit 2
    fi
done

run_orthofinder() {
    local label="$1"
    local input_dir="$2"
    local output_dir="$3"
    local results_dir

    if [ -e "${output_dir}" ]; then
        echo "[ERROR] OrthoFinder output directory already exists for ${label}: ${output_dir}" >&2
        echo "[ERROR] Set RUN_ID to a new value or move the existing directory before resubmitting." >&2
        exit 3
    fi

    orthofinder \
        -f "${input_dir}" \
        -S "${SEARCH_PROGRAM}" \
        -t "${THREADS}" \
        -a "${ANALYSIS_THREADS}" \
        -og \
        -o "${output_dir}"

    results_dir=""
    if [ -s "${output_dir}/Orthogroups/Orthogroups.tsv" ]; then
        results_dir="${output_dir}"
    else
        for candidate in "${output_dir}"/Results* "${output_dir}"/OrthoFinder/Results*; do
            if [ -s "${candidate}/Orthogroups/Orthogroups.tsv" ]; then
                results_dir="${candidate}"
                break
            fi
        done
    fi

    if [ -z "${results_dir}" ] || [ ! -d "${results_dir}" ]; then
        echo "[ERROR] Expected OrthoFinder result directory was not found for ${label} under: ${output_dir}" >&2
        exit 3
    fi

    for required in \
        "${results_dir}/Orthogroups/Orthogroups.tsv" \
        "${results_dir}/Orthogroups/Orthogroups.GeneCount.tsv" \
        "${results_dir}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
    do
        if [ ! -s "${required}" ]; then
            echo "[ERROR] Missing expected OrthoFinder output for ${label}: ${required}" >&2
            exit 3
        fi
    done

    echo "[OK] ${label} OrthoFinder results: ${results_dir}"
    printf '%s\n' "${results_dir}" > "${ANALYSIS_DIR}/orthofinder/results/${label}.latest.txt"
}

run_orthofinder \
    "main" \
    "${ANALYSIS_DIR}/orthofinder/main/pep" \
    "${ANALYSIS_DIR}/orthofinder/results/main_${RUN_ID}"

run_orthofinder \
    "sensitivity" \
    "${ANALYSIS_DIR}/orthofinder/sensitivity/pep" \
    "${ANALYSIS_DIR}/orthofinder/results/sensitivity_${RUN_ID}"
