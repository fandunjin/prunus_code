#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/psal"
ANNOTATION_DIR="${DATA_DIR}/annotation"
GENOME="${ANNOTATION_DIR}/helixer/work/psal.genomic.fa"
HELIXER_GFF="${ANNOTATION_DIR}/helixer/psal.helixer.gff3"
REFINE_DIR="${ANNOTATION_DIR}/eviann"
EVIDENCE_DIR="${REFINE_DIR}/evidence"
POST_DIR="${REFINE_DIR}/postprocess"
PREFIX="psal.eviann"

# Peach is the closest high-quality annotated diploid Prunus reference available.
REFERENCE_LABEL="ppe"
REFERENCE_CDS="${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz"
REFERENCE_PROTEIN="${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz"
EVIDENCE_CDS="${EVIDENCE_DIR}/${REFERENCE_LABEL}.cds.fa"
EVIDENCE_PROTEIN="${EVIDENCE_DIR}/${REFERENCE_LABEL}.protein.fa"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
EVIANN_ENV="${EVIANN_ENV:-eviann}"
PY_ENV="${PY_ENV:-haphic}"
THREADS="${THREADS:-24}"
SWISSPROT="${SWISSPROT:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
POSTPROCESS_SCRIPT="${ROOT_DIR}/scripts/helixer_postprocess_qc.py"
GFF_STATS_SCRIPT="${ROOT_DIR}/scripts/gff3_basic_stats.py"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

copy_fasta() {
    local src="$1"
    local dst="$2"
    if [[ "${src}" == *.gz ]]; then
        gzip -dc "${src}" > "${dst}.tmp"
    else
        cp "${src}" "${dst}.tmp"
    fi
    mv "${dst}.tmp" "${dst}"
}

require_file "${GENOME}" "P. salicina genome prepared by Helixer job"
require_file "${HELIXER_GFF}" "P. salicina Helixer GFF3"
require_file "${REFERENCE_CDS}" "${REFERENCE_LABEL} reference CDS"
require_file "${REFERENCE_PROTEIN}" "${REFERENCE_LABEL} reference protein"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${POSTPROCESS_SCRIPT}" "postprocess/QC script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 statistics script"
mkdir -p "${EVIDENCE_DIR}" "${POST_DIR}"

[ -s "${EVIDENCE_CDS}" ] || copy_fasta "${REFERENCE_CDS}" "${EVIDENCE_CDS}"
[ -s "${EVIDENCE_PROTEIN}" ] || copy_fasta "${REFERENCE_PROTEIN}" "${EVIDENCE_PROTEIN}"
require_file "${EVIDENCE_CDS}" "prepared reference CDS"
require_file "${EVIDENCE_PROTEIN}" "prepared reference protein"

source "${CONDA_SH}"
conda activate "${EVIANN_ENV}"

RAW_GFF="${REFINE_DIR}/psal.genomic.fa.pseudo_label.gff"
FUNCTIONAL_GFF="${REFINE_DIR}/psal.genomic.fa.functional_note.pseudo_label.gff"
if [ ! -s "${RAW_GFF}" ]; then
    cd "${REFINE_DIR}"
    EVIANN_CMD=(
        eviann.sh
        -g "${GENOME}"
        -t "${THREADS}"
        -e "${EVIDENCE_CDS}"
        -p "${EVIDENCE_PROTEIN}"
        -c "${HELIXER_GFF}"
        -f
        --verbose
    )
    if [ -s "${SWISSPROT}" ]; then
        EVIANN_CMD+=(-s "${SWISSPROT}")
    fi
    "${EVIANN_CMD[@]}"
else
    echo "Existing EviAnn GFF found; skipping refinement: ${RAW_GFF}"
fi
require_file "${RAW_GFF}" "EviAnn GFF3"

FINAL_GFF="${RAW_GFF}"
if [ -s "${FUNCTIONAL_GFF}" ]; then
    FINAL_GFF="${FUNCTIONAL_GFF}"
fi

conda activate "${PY_ENV}"
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi
"${PYTHON_BIN}" "${POSTPROCESS_SCRIPT}" \
    --genome "${GENOME}" \
    --gff3 "${FINAL_GFF}" \
    --out-dir "${POST_DIR}" \
    --prefix "${PREFIX}" \
    --short-cds 300 \
    --gene-te-fraction 0.80 \
    --cds-te-fraction 0.50

require_file "${POST_DIR}/${PREFIX}.genes.gff3" "final P. salicina GFF3"
require_file "${POST_DIR}/${PREFIX}.cds.fa" "final P. salicina CDS"
require_file "${POST_DIR}/${PREFIX}.protein.fa" "final P. salicina protein"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${POST_DIR}/${PREFIX}.genes.gff3" --prefix "${PREFIX}" --out-dir "${POST_DIR}/stats"

echo "Final GFF3   : ${POST_DIR}/${PREFIX}.genes.gff3"
echo "Final CDS    : ${POST_DIR}/${PREFIX}.cds.fa"
echo "Final protein: ${POST_DIR}/${PREFIX}.protein.fa"
