#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
#
# Optional overrides:
#   ROOT=/path/to/prunus OUT=/path/to/output THREADS=32 jsub < jobs/psp_extra_qc.sh

set -euo pipefail

ROOT="${ROOT:-${PROJECT_ROOT}/data/psp/psp_raw_qc}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
THREADS="${THREADS:-32}"

[[ -s "${CONDA_SH}" ]] || { echo "ERROR: missing conda init script: ${CONDA_SH}" >&2; exit 1; }
mkdir -p "${OUT}" "${ROOT}/logs"

source "${CONDA_SH}"
conda activate seqkit
command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit unavailable" >&2; exit 1; }

shopt -s nullglob
FILES=(
    "${ROOT}"/data/psp/CRR2018655*.fq.gz
    "${ROOT}"/data/psp/CRR2018655*.fastq.gz
    "${ROOT}"/data/psp/CRR2018656_r2*.fq.gz
    "${ROOT}"/data/psp/CRR2018656_r2*.fastq.gz
)
shopt -u nullglob

(( ${#FILES[@]} > 0 )) || { echo "ERROR: no CRR2018655/CRR2018656_r2 FASTQ files found under ${ROOT}/data/psp" >&2; exit 1; }

echo "=== P. spinulosa extra seqkit statistics ==="
echo "Project root : ${ROOT}"
echo "Input files  : ${#FILES[@]}"
printf '  %s\n' "${FILES[@]}"
echo "Threads      : ${THREADS}"
echo

seqkit stats -a -T -j "${THREADS}" "${FILES[@]}" | tee "${OUT}/psp_extra_seqkit.tsv"

awk -F '\t' '
    NR == 1 {
        print "file\tread_count\ttotal_bases\tread_N50\tQ20_percent\tQ30_percent\tQ1"
        next
    }
    {
        gsub(",", "", $4); gsub(",", "", $5)
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $4, $5, $13, $15, $16, $9
    }
' "${OUT}/psp_extra_seqkit.tsv" > "${OUT}/psp_extra_read_metrics.tsv"

echo
echo "=== extra read metrics ==="
cat "${OUT}/psp_extra_read_metrics.tsv"
echo
echo "Metrics written to: ${OUT}/psp_extra_read_metrics.tsv"
