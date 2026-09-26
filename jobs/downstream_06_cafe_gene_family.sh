#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
CAFE_ENV="${CAFE_ENV:-cafe5}"
THREADS="${THREADS:-16}"
CAFE_DIR="${ANALYSIS_DIR}/cafe5/main"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ORTHOFINDER_LATEST="${ORTHOFINDER_LATEST:-${ANALYSIS_DIR}/orthofinder/results/main.latest.txt}"
SPECIES="${SPECIES:-pzi_A2B2,ppe,pmu_refseq,pdu,pcam,pspe,mdom,msyl,pcom,rch,fve}"
MAX_FAMILY_SIZE="${MAX_FAMILY_SIZE:-100}"
MIN_SPECIES_PRESENT="${MIN_SPECIES_PRESENT:-2}"
LAMBDA_CATEGORIES="${LAMBDA_CATEGORIES:-3}"

mkdir -p "${CAFE_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${CAFE_ENV}"
set -u

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN=python
    else
        echo "[ERROR] Neither ${PYTHON_BIN} nor python is available in CAFE_ENV=${CAFE_ENV}" >&2
        exit 2
    fi
fi

if ! command -v cafe5 >/dev/null 2>&1; then
    echo "[ERROR] Required command not found in CAFE_ENV=${CAFE_ENV}: cafe5" >&2
    exit 2
fi

if [ ! -s "${ORTHOFINDER_LATEST}" ]; then
    echo "[ERROR] Missing OrthoFinder latest pointer: ${ORTHOFINDER_LATEST}" >&2
    echo "[ERROR] Run downstream_03_orthofinder.sh successfully before downstream_06." >&2
    exit 3
fi
read -r of_dir < "${ORTHOFINDER_LATEST}"
if [ ! -d "${of_dir}" ]; then
    echo "[ERROR] OrthoFinder result directory from ${ORTHOFINDER_LATEST} does not exist: ${of_dir}" >&2
    exit 3
fi
if [ ! -s "${of_dir}/Orthogroups/Orthogroups.GeneCount.tsv" ]; then
    echo "[ERROR] Missing OrthoFinder gene-count table: ${of_dir}/Orthogroups/Orthogroups.GeneCount.tsv" >&2
    exit 3
fi
if [ ! -s "${CAFE_DIR}/mcmctree_dated_ultrametric_tree.tre" ]; then
    echo "[ERROR] Missing CAFE dated ultrametric tree: ${CAFE_DIR}/mcmctree_dated_ultrametric_tree.tre" >&2
    echo "[ERROR] Copy the finalized MCMCTree dated tree here before running CAFE5." >&2
    exit 3
fi
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_make_cafe_counts.py" \
    --gene-counts "${of_dir}/Orthogroups/Orthogroups.GeneCount.tsv" \
    --species "${SPECIES}" \
    --out "${CAFE_DIR}/gene_family_counts.tsv" \
    --max-family-size "${MAX_FAMILY_SIZE}" \
    --min-species-present "${MIN_SPECIES_PRESENT}"

cafe5 \
    -i "${CAFE_DIR}/gene_family_counts.tsv" \
    -t "${CAFE_DIR}/mcmctree_dated_ultrametric_tree.tre" \
    -k "${LAMBDA_CATEGORIES}" \
    -c "${THREADS}" \
    -o "${CAFE_DIR}/k${LAMBDA_CATEGORIES}"
