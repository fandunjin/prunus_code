#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
BUSCO_ROOT="${BUSCO_ROOT:-${ROOT_DIR}/haploid_repartition/busco}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

case "${SPECIES}:${LABEL}" in
    pzi:hapA1|pzi:hapA2|pzi:hapB1|pzi:hapB2|psp:hapA1|psp:hapA2|psp:hapB1|psp:hapB2|ppa:hapA|ppa:hapB) ;;
    *) echo "[ERROR] Unsupported SPECIES/LABEL=${SPECIES}/${LABEL}" >&2; exit 2 ;;
esac

PROTEIN="${PROTEIN:-${INPUT_ROOT}/${SPECIES}/pep/${SPECIES}_${LABEL}.pep.fa}"
RUN_NAME="${SPECIES}_${LABEL}_eviann_protein_eudicots"
mkdir -p "${BUSCO_ROOT}" "${ROOT_DIR}/logs"
[[ -s "${PROTEIN}" ]] || { echo "[ERROR] Missing EviAnn protein: ${PROTEIN}" >&2; exit 2; }
[[ -s "${BUSCO_LINEAGE}" ]] || { echo "[ERROR] Missing BUSCO lineage: ${BUSCO_LINEAGE}" >&2; exit 2; }
[[ -s "${CONDA_SH}" ]] || { echo "[ERROR] Missing conda initialization: ${CONDA_SH}" >&2; exit 2; }

source "${CONDA_SH}"
conda activate "${BUSCO_ENV}"
cd "${BUSCO_ROOT}"
busco -i "${PROTEIN}" -m proteins -l "${BUSCO_LINEAGE}" -c "${THREADS}" --offline -o "${RUN_NAME}" --force
find "${BUSCO_ROOT}/${RUN_NAME}" -name 'short_summary*.txt' -type f -size +0 -print | sort
echo "[OK] BUSCO output: ${BUSCO_ROOT}/${RUN_NAME}"
