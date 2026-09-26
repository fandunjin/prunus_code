#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/psp/haploid/genome/psp_A1A2B1B2.renamed.chr32.fa}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/data/psp/annotation/busco_genome}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
RUN_NAME="${RUN_NAME:-psp32_current_genome_eudicotyledons_odb12}"
FORCE="${FORCE:-0}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
case "${FORCE}" in
    0|1) ;;
    *) die "FORCE must be 0 or 1 (observed: ${FORCE})" ;;
esac

require_file "${PSP_GENOME}" "PSP 32-chromosome genome FASTA"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

set +u
source "${CONDA_SH}"
set -u
conda activate "${BUSCO_ENV}"
require_cmd busco

RUN_DIR="${OUT_DIR}/${RUN_NAME}"
SUMMARY="$(find "${RUN_DIR}" -maxdepth 3 -type f -name 'short_summary*.txt' -size +0c 2>/dev/null | head -n 1 || true)"
if [[ -n "${SUMMARY}" && "${FORCE}" != "1" ]]; then
    echo "Existing complete BUSCO result found; skipping: ${SUMMARY}"
    exit 0
fi

echo "=== BUSCO: current PSP 32-chromosome genome ==="
echo "Genome  : ${PSP_GENOME}"
echo "Lineage : ${BUSCO_LINEAGE}"
echo "Mode    : genome"
echo "Threads : ${THREADS}"
echo "Output  : ${RUN_DIR}"
echo "Force   : ${FORCE}"

BUSCO_ARGS=(
    -i "${PSP_GENOME}"
    -m genome
    -l "${BUSCO_LINEAGE}"
    -c "${THREADS}"
    --offline
    -o "${RUN_NAME}"
    --out_path "${OUT_DIR}"
)
if [[ "${FORCE}" == "1" ]]; then
    BUSCO_ARGS+=(--force)
elif [[ -d "${RUN_DIR}" ]]; then
    BUSCO_ARGS+=(--restart)
fi

busco "${BUSCO_ARGS[@]}"

SUMMARY="$(find "${RUN_DIR}" -maxdepth 3 -type f -name 'short_summary*.txt' -size +0c | head -n 1 || true)"
[[ -n "${SUMMARY}" ]] || die "BUSCO completed without a non-empty short summary in ${RUN_DIR}"
echo "BUSCO completed. Summary: ${SUMMARY}"
