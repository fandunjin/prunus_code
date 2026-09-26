#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# BUSCO jobs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/ppa}"
PPA_EVIANN_DIR="${PPA_EVIANN_DIR:-${DATA_DIR}/annotation/eviann}"
REFERENCE_MODE="${REFERENCE_MODE:-psp_pzi}"
if [[ "${REFERENCE_MODE}" == "ppe" ]]; then
    DEFAULT_POST_DIR="${PPA_EVIANN_DIR}/postprocess_ppe"
    DEFAULT_PREFIX="ppa.ppe.eviann"
    DEFAULT_RUN_NAME="ppa_ppe_eviann_protein_eudicots"
elif [[ "${REFERENCE_MODE}" == "legacy" ]]; then
    DEFAULT_POST_DIR="${PPA_EVIANN_DIR}/postprocess_legacy"
    DEFAULT_PREFIX="ppa.legacy.eviann"
    DEFAULT_RUN_NAME="ppa_legacy_eviann_protein_eudicots"
else
    DEFAULT_POST_DIR="${PPA_EVIANN_DIR}/postprocess"
    DEFAULT_PREFIX="ppa.eviann"
    DEFAULT_RUN_NAME="ppa_eviann_protein_eudicots"
fi
PROTEIN="${PROTEIN:-${DEFAULT_POST_DIR}/${DEFAULT_PREFIX}.protein.fa}"
GFF3="${GFF3:-${DEFAULT_POST_DIR}/${DEFAULT_PREFIX}.genes.gff3}"
OUT_BASE="${OUT_BASE:-${PPA_EVIANN_DIR}/busco}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
RUN_NAME="${RUN_NAME:-${DEFAULT_RUN_NAME}}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

case "${REFERENCE_MODE}" in
    legacy|psp_pzi|ppe) ;;
    *) die "REFERENCE_MODE must be legacy, psp_pzi or ppe (observed: ${REFERENCE_MODE})" ;;
esac

# Older PPA runs wrote their complete annotation to data/ppa/postprocess rather
# than data/ppa/annotation/eviann/postprocess.  Accept that legacy result only
# when the current output is absent, and keep the selected source in the log.
LEGACY_PROTEIN="${DATA_DIR}/postprocess/ppa.eviann.protein.fa"
LEGACY_GFF3="${DATA_DIR}/postprocess/ppa.eviann.genes.gff3"
if [[ "${REFERENCE_MODE}" == "legacy" && ! -s "${PROTEIN}" && -s "${LEGACY_PROTEIN}" ]]; then
    echo "WARNING: current PPA EviAnn protein path is absent; using legacy output: ${LEGACY_PROTEIN}" >&2
    PROTEIN="${LEGACY_PROTEIN}"
    GFF3="${LEGACY_GFF3}"
    OUT_BASE="${DATA_DIR}/postprocess/busco"
fi
require_file "${PROTEIN}" "EviAnn final PPA protein FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
mkdir -p "${OUT_BASE}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${BUSCO_ENV}"
require_cmd busco

cd "${OUT_BASE}"
echo "=== BUSCO: final PPA EviAnn protein annotation ==="
echo "Reference mode: ${REFERENCE_MODE}"
echo "Protein : ${PROTEIN}"
if [[ -s "${GFF3}" ]]; then
    echo "GFF3    : ${GFF3}"
else
    echo "GFF3    : not found (protein-only validation)"
fi
PROTEIN_RECORDS="$(awk '/^>/{n++} END{print n+0}' "${PROTEIN}")"
echo "Protein records: ${PROTEIN_RECORDS}"
if [[ -s "${GFF3}" ]]; then
    GFF_GENES="$(awk -F '\t' '!/^#/ && $3 == "gene" {n++} END{print n+0}' "${GFF3}")"
    GFF_MRNA="$(awk -F '\t' '!/^#/ && ($3 == "mRNA" || $3 == "transcript") {n++} END{print n+0}' "${GFF3}")"
    GFF_CDS="$(awk -F '\t' '!/^#/ && $3 == "CDS" {n++} END{print n+0}' "${GFF3}")"
    echo "GFF3 records: genes=${GFF_GENES}, transcripts=${GFF_MRNA}, CDS_parts=${GFF_CDS}"
fi
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
