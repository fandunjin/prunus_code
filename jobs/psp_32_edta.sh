#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
INPUT="${INPUT:-${DATA_DIR}/haphic/merge/04.build/final_assembly_qc_outputs/psp_2n32.chromosomes.fa}"
OUT_DIR="${OUT_DIR:-${DATA_DIR}/annotation/repeats/chr32_edta}"
WORK_FASTA="${WORK_FASTA:-${OUT_DIR}/psp_2n32.chromosomes.fa}"
THREADS="${THREADS:-24}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${INPUT}" "PSP 2n=32 chromosome FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
if [ ! -s "${EDTA_CDS}" ] && [ -s "${EDTA_CDS}.gz" ]; then
    EDTA_CDS="${EDTA_CDS}.gz"
fi

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"
if [ ! -s "${WORK_FASTA}" ]; then
    cp -p "${INPUT}" "${WORK_FASTA}"
fi

source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v EDTA.pl >/dev/null 2>&1 || {
    echo "ERROR: EDTA.pl not found in environment ${EDTA_ENV}" >&2
    exit 1
}

cd "${OUT_DIR}"
echo "=== EDTA repeat annotation: PSP 2n=32 chromosomes ==="
echo "Input : ${WORK_FASTA}"
echo "Out   : ${OUT_DIR}"
echo "Force : ${EDTA_FORCE}"
echo "Overwrite existing EDTA intermediates: ${EDTA_OVERWRITE}"

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

require_file "${WORK_FASTA}.mod.EDTA.TEanno.gff3" "PSP EDTA TE annotation GFF3"
require_file "${WORK_FASTA}.mod.EDTA.TElib.fa" "PSP EDTA TE library"
require_file "${WORK_FASTA}.mod.MAKER.masked" "PSP EDTA masked genome"
echo "PSP chr32 EDTA completed: ${OUT_DIR}"
