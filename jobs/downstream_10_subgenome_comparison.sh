#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
HAP_GENOME_DIR="${HAP_RENAME_DIR}/genome"
HAP_CDS_DIR="${HAP_RENAME_DIR}/cds"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
THREADS="${THREADS:-32}"
FORCE="${FORCE:-0}"
SUB_DIR="${ANALYSIS_DIR}/subgenomes"

mkdir -p "${SUB_DIR}" "${ANALYSIS_DIR}/logs"
[[ -s "${PACKAGE_CHECKSUMS}" ]] || { echo "[ERROR] Missing renamed-package checksums: ${PACKAGE_CHECKSUMS}" >&2; exit 2; }

set +u
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
set -u

for exe in EDTA.pl LTR_retriever RepeatModeler RepeatMasker gt; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in EDTA_ENV=${EDTA_ENV}: ${exe}" >&2
        exit 2
    fi
done

for label in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2; do
    genome="${HAP_GENOME_DIR}/${label}.fa"
    cds="${HAP_CDS_DIR}/${label}.cds.fa"
    [[ -s "${genome}" ]] || { echo "[ERROR] Missing renamed genome: ${genome}" >&2; exit 2; }
    [[ -s "${cds}" ]] || { echo "[ERROR] Missing renamed CDS: ${cds}" >&2; exit 2; }
    run_dir="${SUB_DIR}/edta/${label}"
    saved_checksums="${run_dir}/.haploid_rename.checksums.md5"
    overwrite="${FORCE}"
    mkdir -p "${run_dir}"
    if ! cmp -s "${PACKAGE_CHECKSUMS}" "${saved_checksums}"; then
        overwrite=1
        echo "[INFO] ${label}: renamed package changed; rebuilding EDTA outputs."
    fi
    cd "${run_dir}"
    EDTA.pl \
        --genome "${genome}" \
        --cds "${cds}" \
        --species others \
        --step all \
        --sensitive 1 \
        --anno 1 \
        --threads "${THREADS}" \
        --overwrite "${overwrite}"
    cp -f "${PACKAGE_CHECKSUMS}" "${saved_checksums}"
done

echo "Combine EDTA summaries, LTR age tables, BUSCO, gene density, and homeolog Ks/4DTV downstream."
