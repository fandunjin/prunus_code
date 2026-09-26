#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Set PPA_GENOME to evaluate another coordinate-compatible PPA FASTA instead.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/ppa}"
PPA_GENOME="${PPA_GENOME:-${DATA_DIR}/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
OUT_DIR="${OUT_DIR:-${DATA_DIR}/busco_genome}"
INPUT_DIR="${INPUT_DIR:-${OUT_DIR}/input}"
INPUT_FASTA="${INPUT_FASTA:-${INPUT_DIR}/ppa.genome.fa}"
RUN_NAME="${RUN_NAME:-ppa_genome_eudicotyledons_odb12}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
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

case "${FORCE}" in
    0|1) ;;
    *) die "FORCE must be 0 or 1 (observed: ${FORCE})" ;;
esac
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"

require_file "${PPA_GENOME}" "PPA genome FASTA (or gzip-compressed FASTA)"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${INPUT_DIR}" "${ROOT_DIR}/logs"

# BUSCO accepts FASTA input, but staging a plain FASTA makes the run
# independent of version-specific handling of gzip-compressed input.
if [[ "${PPA_GENOME}" == *.gz ]]; then
    require_cmd gzip
    gzip -t "${PPA_GENOME}" || die "Invalid gzip file: ${PPA_GENOME}"
    if [[ ! -s "${INPUT_FASTA}" || "${FORCE}" == "1" ]]; then
        tmp="${INPUT_FASTA}.tmp.$$"
        rm -f "${tmp}"
        gzip -dc "${PPA_GENOME}" > "${tmp}"
        mv -f "${tmp}" "${INPUT_FASTA}"
    fi
else
    INPUT_FASTA="${PPA_GENOME}"
fi
require_file "${INPUT_FASTA}" "staged PPA genome FASTA"

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

echo "=== BUSCO: PPA genome mode ==="
echo "Source genome : ${PPA_GENOME}"
echo "Input FASTA   : ${INPUT_FASTA}"
echo "Lineage       : ${BUSCO_LINEAGE}"
echo "Mode          : genome"
echo "Threads       : ${THREADS}"
echo "Output        : ${RUN_DIR}"
echo "Force         : ${FORCE}"

BUSCO_ARGS=(
    -i "${INPUT_FASTA}"
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
