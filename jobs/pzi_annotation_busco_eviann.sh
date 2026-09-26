#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_BASE="${OUT_BASE:-${ROOT_DIR}/annotation/genes/pzi_eviann/busco}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
RUN_NAME="${RUN_NAME:-pzi32_eviann_protein_eudicots}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

require_file "${PROTEIN}" "EviAnn final protein FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
mkdir -p "${OUT_BASE}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${BUSCO_ENV}"
require_cmd busco

cd "${OUT_BASE}"
echo "=== BUSCO: final PZI EviAnn protein annotation ==="
echo "Protein : ${PROTEIN}"
echo "Lineage : ${BUSCO_LINEAGE}"
echo "Output  : ${OUT_BASE}/${RUN_NAME}"
echo "Threads : ${THREADS}"

busco \
    -i "${PROTEIN}" \
    -m proteins \
    -l "${BUSCO_LINEAGE}" \
    -c "${THREADS}" \
    --offline \
    -o "${RUN_NAME}" \
    --force

echo "BUSCO completed. Summary files:"
find "${OUT_BASE}/${RUN_NAME}" -name "short_summary*.txt" -type f -print
