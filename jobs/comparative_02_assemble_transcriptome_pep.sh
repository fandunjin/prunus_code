#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/transcriptome_accessions.tsv"
TX_DIR="${ROOT_DIR}/data/transcriptomes"
OUT_DIR="${ROOT_DIR}/comparative_genomics/transcriptome_derived_pep"
THREADS=32
MEMORY_GB=180

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

if [[ ! -s "${ACCESSIONS}" ]]; then
    echo "ERROR: missing ${ACCESSIONS}"
    exit 1
fi

if command -v module >/dev/null 2>&1; then
    module load trinity 2>/dev/null || true
    module load transdecoder 2>/dev/null || true
fi

if ! command -v Trinity >/dev/null 2>&1; then
    echo "ERROR: Trinity not found. Load Trinity or activate an environment containing Trinity."
    exit 1
fi
if ! command -v TransDecoder.LongOrfs >/dev/null 2>&1; then
    echo "ERROR: TransDecoder.LongOrfs not found."
    exit 1
fi
if ! command -v TransDecoder.Predict >/dev/null 2>&1; then
    echo "ERROR: TransDecoder.Predict not found."
    exit 1
fi

tail -n +2 "${ACCESSIONS}" | while IFS=$'\t' read -r label species source accession read_type tissue note; do
    [[ -z "${label}" || -z "${accession}" ]] && continue
    [[ "${source}" != "SRA" ]] && continue

    fastq_dir="${TX_DIR}/${label}/fastq"
    work_dir="${OUT_DIR}/${label}"
    mkdir -p "${work_dir}"

    left=$(find "${fastq_dir}" -maxdepth 1 -name "*_1.fastq.gz" | sort | paste -sd, -)
    right=$(find "${fastq_dir}" -maxdepth 1 -name "*_2.fastq.gz" | sort | paste -sd, -)
    single=$(find "${fastq_dir}" -maxdepth 1 -name "*.fastq.gz" ! -name "*_1.fastq.gz" ! -name "*_2.fastq.gz" | sort | paste -sd, -)

    echo "=== ${label} ${species} ==="
    if [[ -n "${left}" && -n "${right}" ]]; then
        Trinity --seqType fq --left "${left}" --right "${right}" \
            --CPU "${THREADS}" --max_memory "${MEMORY_GB}G" \
            --output "${work_dir}/trinity"
    elif [[ -n "${single}" ]]; then
        Trinity --seqType fq --single "${single}" \
            --CPU "${THREADS}" --max_memory "${MEMORY_GB}G" \
            --output "${work_dir}/trinity"
    else
        echo "ERROR: no FASTQ files found for ${label} in ${fastq_dir}"
        exit 1
    fi

    transcript_fa="${work_dir}/trinity/Trinity.fasta"
    TransDecoder.LongOrfs -t "${transcript_fa}"
    TransDecoder.Predict -t "${transcript_fa}"

    cp "${transcript_fa}.transdecoder.pep" "${OUT_DIR}/${label}.transdecoder.pep"
    cp "${transcript_fa}.transdecoder.cds" "${OUT_DIR}/${label}.transdecoder.cds"
done

echo "Transcriptome-derived proteins written to ${OUT_DIR}"
