#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${ROOT_DIR}/annotation/chr32_edta"
WORK_FASTA="${OUT_DIR}/pzi_2n32.chromosomes.fa"
THREADS=24
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
EDTA_ENV="EDTA"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${INPUT_FASTA}" "2n=32 chromosome-only FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"
if [ ! -s "${WORK_FASTA}" ]; then
    cp -p "${INPUT_FASTA}" "${WORK_FASTA}"
fi

source "${CONDA_SH}"
conda activate "${EDTA_ENV}"

cd "${OUT_DIR}"

echo "=== EDTA repeat annotation: 2n=32 chromosomes ==="
echo "Input : ${WORK_FASTA}"
echo "Out   : ${OUT_DIR}"
echo "Force : ${EDTA_FORCE}"
echo "Overwrite existing EDTA intermediates: ${EDTA_OVERWRITE}"
echo

if [ -s "${EDTA_CDS}" ]; then
    echo "CDS filter: ${EDTA_CDS}"
else
    echo "CDS filter: not used; CDS file was not found."
fi

if [ ! -s "${WORK_FASTA}.mod.EDTA.TEanno.gff3" ]; then
    EDTA_CMD=(
        EDTA.pl
        --genome "${WORK_FASTA}"
        --species others
        --step all
        --anno 1
        --sensitive 1
        --force "${EDTA_FORCE}"
        --overwrite "${EDTA_OVERWRITE}"
        --threads "${THREADS}"
    )
    if [ -s "${EDTA_CDS}" ]; then
        EDTA_CMD+=(--cds "${EDTA_CDS}")
    fi
    "${EDTA_CMD[@]}"
else
    echo "Existing EDTA TE annotation found, skipping EDTA: ${WORK_FASTA}.mod.EDTA.TEanno.gff3"
fi

require_file "${WORK_FASTA}.mod.EDTA.TEanno.gff3" "EDTA TE annotation GFF3"
require_file "${WORK_FASTA}.mod.EDTA.TElib.fa" "EDTA TE library"

echo "EDTA chr32 completed."
