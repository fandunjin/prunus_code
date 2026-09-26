#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/transcriptome_accessions.tsv"
OUT_DIR="${ROOT_DIR}/data/transcriptomes"

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

if [[ ! -s "${ACCESSIONS}" ]]; then
    echo "ERROR: missing ${ACCESSIONS}"
    echo "Copy data/transcriptome_accessions_template.tsv to data/transcriptome_accessions.tsv and fill accessions first."
    exit 1
fi

if command -v module >/dev/null 2>&1; then
    module load sra-tools 2>/dev/null || true
fi

if ! command -v prefetch >/dev/null 2>&1; then
    echo "ERROR: prefetch not found. Load sra-tools first."
    exit 1
fi
if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "ERROR: fasterq-dump not found. Load sra-tools first."
    exit 1
fi

tail -n +2 "${ACCESSIONS}" | while IFS=$'\t' read -r label species source accession read_type tissue note; do
    [[ -z "${label}" || -z "${accession}" ]] && continue
    if [[ "${source}" != "SRA" ]]; then
        echo "SKIP ${label}: source=${source}, only SRA is handled by this script."
        continue
    fi

    sample_dir="${OUT_DIR}/${label}"
    mkdir -p "${sample_dir}/sra" "${sample_dir}/fastq"
    echo "=== ${label} ${species} ${accession} ==="
    prefetch "${accession}" --output-directory "${sample_dir}/sra"
    fasterq-dump "${sample_dir}/sra/${accession}" \
        --outdir "${sample_dir}/fastq" \
        --threads 8 \
        --split-files
    gzip -f "${sample_dir}/fastq/"*.fastq
done

echo "Transcriptome FASTQ files written to ${OUT_DIR}"
