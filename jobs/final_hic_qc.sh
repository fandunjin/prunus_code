#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
CHROM_SIZES="${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.sizes"

HIC_R1="${PROJECT_ROOT}/HIC/H260100013A_1.fq.gz"
HIC_R2="${PROJECT_ROOT}/HIC/H260100013A_2.fq.gz"

OUT_DIR="${ROOT_DIR}/final_assembly/hic_qc"

THREADS=48
SAMTOOLS_THREADS=12
FILTER_THREADS=36

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
HAPHIC_ENV="haphic"
HAPHIC_DIR="${HOME}/cx/gitsoftware/HapHiC"
CONTACT_MATRIX_SCRIPT="${ROOT_DIR}/scripts/hic_contact_matrix.py"

# Minimum mapping quality used only for the chromosome-pair summary.
MIN_MAPQ_FOR_CONTACTS=10

# ===================== Pipeline starts here =====================

BAM="${OUT_DIR}/HiC.final.bam"
FILTERED_BAM="${OUT_DIR}/HiC.final.filtered.bam"
SORTED_BAM="${OUT_DIR}/HiC.final.filtered.sorted.bam"
CONTACT_TSV="${OUT_DIR}/chr_pair_contacts.tsv"
MATRIX_TSV="${OUT_DIR}/chr_pair_contact_matrix.tsv"
SUMMARY_TSV="${OUT_DIR}/chr_contact_support.tsv"
HEATMAP_PNG="${OUT_DIR}/chr_pair_contact_heatmap.png"

require_file() {
    local file_path="$1"
    local label="$2"

    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${FINAL_FASTA}" "Final FASTA"
require_file "${CHROM_SIZES}" "Chromosome sizes"
require_file "${HIC_R1}" "Hi-C R1"
require_file "${HIC_R2}" "Hi-C R2"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${CONTACT_MATRIX_SCRIPT}" "Hi-C contact matrix Python script"

mkdir -p "${OUT_DIR}"

source "${CONDA_SH}"
conda activate "${HAPHIC_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo "=== Final Hi-C support QC ==="
echo "Final FASTA : ${FINAL_FASTA}"
echo "Hi-C R1     : ${HIC_R1}"
echo "Hi-C R2     : ${HIC_R2}"
echo "Output dir  : ${OUT_DIR}"
echo

if [ ! -s "${FINAL_FASTA}.bwt" ]; then
    echo ">>> Indexing final FASTA with bwa"
    bwa index "${FINAL_FASTA}"
else
    echo ">>> bwa index exists, skipping"
fi

if [ ! -s "${BAM}" ]; then
    echo ">>> Aligning Hi-C reads to final FASTA"
    bwa mem -5SP -t "${THREADS}" "${FINAL_FASTA}" "${HIC_R1}" "${HIC_R2}" \
        | samblaster \
        | samtools view - -@ "${SAMTOOLS_THREADS}" -S -h -b -F 3340 -o "${BAM}"
else
    echo ">>> ${BAM} exists, skipping alignment"
fi

samtools flagstat -@ "${SAMTOOLS_THREADS}" "${BAM}" > "${OUT_DIR}/HiC.final.flagstat.txt"
samtools stats -@ "${SAMTOOLS_THREADS}" "${BAM}" > "${OUT_DIR}/HiC.final.samtools_stats.txt"

if [ ! -s "${FILTERED_BAM}" ]; then
    echo ">>> Filtering Hi-C BAM with HapHiC filter_bam"
    "${HAPHIC_DIR}/utils/filter_bam" "${BAM}" 1 --nm 3 --threads "${FILTER_THREADS}" \
        | samtools view - -b -@ "${SAMTOOLS_THREADS}" -o "${FILTERED_BAM}"
else
    echo ">>> ${FILTERED_BAM} exists, skipping HapHiC filtering"
fi

if [ ! -s "${SORTED_BAM}" ]; then
    echo ">>> Sorting and indexing filtered BAM"
    samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" "${FILTERED_BAM}"
    samtools index -@ "${THREADS}" "${SORTED_BAM}"
else
    echo ">>> ${SORTED_BAM} exists, skipping sort"
fi

samtools flagstat -@ "${SAMTOOLS_THREADS}" "${FILTERED_BAM}" > "${OUT_DIR}/HiC.final.filtered.flagstat.txt"
samtools idxstats "${SORTED_BAM}" > "${OUT_DIR}/HiC.final.filtered.idxstats.txt"

echo ">>> Building chromosome-pair contact matrix"
"${PYTHON_BIN}" "${CONTACT_MATRIX_SCRIPT}" \
    --bam "${SORTED_BAM}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --contact-tsv "${CONTACT_TSV}" \
    --matrix-tsv "${MATRIX_TSV}" \
    --summary-tsv "${SUMMARY_TSV}" \
    --heatmap-png "${HEATMAP_PNG}" \
    --min-mapq "${MIN_MAPQ_FOR_CONTACTS}"

{
    printf "Final Hi-C QC manual checks\n"
    printf "===========================\n\n"
    printf "1. Inspect %s or %s.\n" "${HEATMAP_PNG}" "${MATRIX_TSV}"
    printf "2. Each chromosome should have strong cis signal and weaker off-diagonal signal.\n"
    printf "3. Review %s; chromosomes with unusually low cis_percent need manual contact-map inspection.\n" "${SUMMARY_TSV}"
    printf "4. Keep Juicebox before/after screenshots beside this directory for the final assembly record.\n"
} > "${OUT_DIR}/manual_review_checklist.txt"

echo
echo "=== Final Hi-C QC outputs ==="
echo "Raw BAM flagstat       : ${OUT_DIR}/HiC.final.flagstat.txt"
echo "Filtered BAM flagstat  : ${OUT_DIR}/HiC.final.filtered.flagstat.txt"
echo "Filtered BAM idxstats  : ${OUT_DIR}/HiC.final.filtered.idxstats.txt"
echo "Contact TSV           : ${CONTACT_TSV}"
echo "Contact matrix        : ${MATRIX_TSV}"
echo "Chromosome support    : ${SUMMARY_TSV}"
echo "Heatmap               : ${HEATMAP_PNG}"
