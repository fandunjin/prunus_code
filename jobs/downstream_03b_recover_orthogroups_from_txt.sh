#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PY_ENV="${PY_ENV:-haphic}"
PYTHON_BIN="${PYTHON_BIN:-python}"
RUN_DIR="${RUN_DIR:-${ANALYSIS_DIR}/orthofinder/results/main_387991/Results_Jun28}"
SET_NAME="${SET_NAME:-main}"

mkdir -p "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${PY_ENV}"
set -u

orthogroups_txt="${RUN_DIR}/Orthogroups/Orthogroups.txt"
pep_manifest="${ANALYSIS_DIR}/orthofinder/${SET_NAME}/pep_manifest.tsv"
working_dir="${RUN_DIR}/WorkingDirectory"
out_dir="${RUN_DIR}/Orthogroups"

if [ ! -s "${orthogroups_txt}" ]; then
    echo "[ERROR] Missing Orthogroups.txt: ${orthogroups_txt}" >&2
    exit 3
fi
if [ ! -s "${pep_manifest}" ]; then
    echo "[ERROR] Missing protein manifest: ${pep_manifest}" >&2
    exit 3
fi

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_recover_orthogroups_tables.py" \
    --orthogroups-txt "${orthogroups_txt}" \
    --pep-manifest "${pep_manifest}" \
    --working-directory "${working_dir}" \
    --out-dir "${out_dir}"

for required in \
    "${out_dir}/Orthogroups.tsv" \
    "${out_dir}/Orthogroups.GeneCount.tsv" \
    "${out_dir}/Orthogroups_SingleCopyOrthologues.txt"
do
    if [ ! -s "${required}" ]; then
        echo "[ERROR] Recovery did not create expected output: ${required}" >&2
        exit 3
    fi
done

printf '%s\n' "${RUN_DIR}" > "${ANALYSIS_DIR}/orthofinder/results/${SET_NAME}.latest.txt"
echo "[OK] Recovered OrthoFinder tables for ${SET_NAME}: ${RUN_DIR}"
echo "[OK] Wrote pointer: ${ANALYSIS_DIR}/orthofinder/results/${SET_NAME}.latest.txt"
