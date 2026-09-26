#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/comparative_genomics_genome_download_list.tsv"
OUT_DIR="${ROOT_DIR}/comparative_genomics/genomes/raw_downloads"
REPORT="${OUT_DIR}/download_report.tsv"

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"
cd "${ROOT_DIR}"

if ! command -v datasets >/dev/null 2>&1; then
  echo "ERROR: NCBI datasets CLI not found in PATH." >&2
  echo "Install or load datasets first, then rerun this job." >&2
  exit 1
fi

echo -e "label\tspecies\tquery\tstatus\tzip_path" > "${REPORT}"

tail -n +2 "${LIST}" | awk -F '\t' '$1=="core" || $1=="outgroup" {print}' | while IFS=$'\t' read -r set priority label species query role use_in_main avoid_reason notes; do
  if [ "${label}" = "pzi_copy1" ]; then
    echo -e "${label}\t${species}\tlocal\tSKIP_LOCAL\tfinal_genome/pzi_copy1.fa" >> "${REPORT}"
    continue
  fi

  ZIP="${OUT_DIR}/${label}.ncbi_dataset.zip"
  echo "=== Downloading ${label}: ${species} (${query}) ==="

  if [[ "${query}" == GCA_* || "${query}" == GCF_* ]]; then
    datasets download genome accession "${query}" \
      --include genome,gff3,cds,protein,seq-report \
      --filename "${ZIP}"
  else
    datasets download genome taxon "${query}" \
      --reference \
      --assembly-version current \
      --include genome,gff3,cds,protein,seq-report \
      --filename "${ZIP}" || \
    datasets download genome taxon "${query}" \
      --assembly-level chromosome,complete \
      --assembly-version current \
      --include genome,gff3,cds,protein,seq-report \
      --filename "${ZIP}"
  fi

  if [ -s "${ZIP}" ]; then
    echo -e "${label}\t${species}\t${query}\tDOWNLOADED\t${ZIP}" >> "${REPORT}"
  else
    echo -e "${label}\t${species}\t${query}\tFAILED_OR_EMPTY\t${ZIP}" >> "${REPORT}"
  fi
done

echo "Download report: ${REPORT}"
