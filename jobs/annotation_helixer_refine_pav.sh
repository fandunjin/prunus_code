#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/pav"
ANNOTATION_DIR="${DATA_DIR}/annotation"
GENOME="${ANNOTATION_DIR}/helixer/work/pav.genomic.fa"
HELIXER_GFF="${ANNOTATION_DIR}/helixer/pav.helixer.gff3"
REFINE_DIR="${ANNOTATION_DIR}/eviann"
EVIDENCE_DIR="${REFINE_DIR}/evidence"
POST_DIR="${REFINE_DIR}/postprocess"
PREFIX="pav.eviann"

# Same-species RefSeq annotation is evidence only; coordinates remain GCA_014155035.1.
REFERENCE_LABEL="pav_GCF002207925"
REFERENCE_CDS="${REFERENCE_CDS:-}"
REFERENCE_PROTEIN="${REFERENCE_PROTEIN:-}"
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

find_reference_file() {
    local pattern="$1"
    find "${ROOT_DIR}/data" -type f -name "${pattern}" -print -quit
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

if [ -z "${REFERENCE_CDS}" ]; then
    REFERENCE_CDS="$(find_reference_file 'GCF_002207925.1_PAV_r1.0_cds_from_genomic.fna.gz')"
fi
if [ -z "${REFERENCE_PROTEIN}" ]; then
    REFERENCE_PROTEIN="$(find_reference_file 'GCF_002207925.1_PAV_r1.0_protein.faa.gz')"
fi

require_file "${GENOME}" "P. avium GCA_014155035.1 genome prepared by Helixer job"
require_file "${HELIXER_GFF}" "P. avium GCA_014155035.1 Helixer GFF3"
require_file "${REFERENCE_CDS}" "GCF_002207925.1 reference CDS"
require_file "${REFERENCE_PROTEIN}" "GCF_002207925.1 reference protein"
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

RAW_GFF="${REFINE_DIR}/pav.genomic.fa.pseudo_label.gff"
FUNCTIONAL_GFF="${REFINE_DIR}/pav.genomic.fa.functional_note.pseudo_label.gff"
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

require_file "${POST_DIR}/${PREFIX}.genes.gff3" "final P. avium GFF3"
require_file "${POST_DIR}/${PREFIX}.cds.fa" "final P. avium CDS"
require_file "${POST_DIR}/${PREFIX}.protein.fa" "final P. avium protein"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${POST_DIR}/${PREFIX}.genes.gff3" --prefix "${PREFIX}" --out-dir "${POST_DIR}/stats"

echo "Final GFF3   : ${POST_DIR}/${PREFIX}.genes.gff3"
echo "Final CDS    : ${POST_DIR}/${PREFIX}.cds.fa"
echo "Final protein: ${POST_DIR}/${PREFIX}.protein.fa"
