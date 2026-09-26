#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
COPY2_FASTA="${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy2.hap_provenance.chr16.fa"
WORK_DIR="${ROOT_DIR}/final_assembly/copy1_copy2_syri"
THREADS=24
MINIMAP_PRESET="asm5"
CHR_COUNT=16

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SYRI_ENV="syri_run"
HELPER_SCRIPT="${ROOT_DIR}/scripts/haplotype_syri_helpers.py"

SELF_PAF="${WORK_DIR}/pzi_copy1_vs_copy2.raw.paf"
SELF_SUMMARY="${WORK_DIR}/pzi_copy1_vs_copy2.raw.same_chr_alignment.tsv"
COPY2_ORIENTED="${WORK_DIR}/pzi_copy2.oriented_to_copy1.chr16.fa"
ORIENTATION_TSV="${WORK_DIR}/pzi_copy2.oriented_to_copy1.orientation.tsv"
ALN_SAM="${WORK_DIR}/pzi_copy1_vs_copy2.oriented.sam"
ALN_BAM="${WORK_DIR}/pzi_copy1_vs_copy2.oriented.bam"
SYRI_PREFIX="pzi_copy1_vs_copy2.oriented."
SYRI_OUT="${WORK_DIR}/${SYRI_PREFIX}syri.out"
GENOMES_TXT="${WORK_DIR}/pzi_copy1_vs_copy2.genomes.txt"

FORCE_RERUN="${FORCE_RERUN:-0}"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

remove_if_forced() {
    if [ "${FORCE_RERUN}" -eq 1 ]; then
        rm -f "$@"
    fi
}

require_file "${COPY1_FASTA}" "copy1 n=16 FASTA"
require_file "${COPY2_FASTA}" "copy2 n=16 FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${HELPER_SCRIPT}" "haplotype SyRI helper script"

mkdir -p "${WORK_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${SYRI_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== copy1 vs copy2 haplotype synteny/SyRI ==="
echo "copy1 reference: ${COPY1_FASTA}"
echo "copy2 query    : ${COPY2_FASTA}"
echo "work dir       : ${WORK_DIR}"
echo "threads        : ${THREADS}"
echo

cd "${WORK_DIR}"

remove_if_forced "${SELF_PAF}" "${SELF_SUMMARY}" "${COPY2_ORIENTED}" "${ORIENTATION_TSV}"
if [ ! -s "${SELF_PAF}" ]; then
    echo ">>> Minimap2 PAF: raw copy2 against copy1"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${COPY1_FASTA}" "${COPY2_FASTA}" > "${SELF_PAF}"
else
    echo ">>> Existing raw PAF found, skipping: ${SELF_PAF}"
fi

"${PYTHON_BIN}" "${HELPER_SCRIPT}" summarize-paf \
    --paf "${SELF_PAF}" \
    --out-tsv "${SELF_SUMMARY}" \
    --chr-count "${CHR_COUNT}"

"${PYTHON_BIN}" "${HELPER_SCRIPT}" orient-query \
    --paf "${SELF_PAF}" \
    --query-in "${COPY2_FASTA}" \
    --query-out "${COPY2_ORIENTED}" \
    --orientation-tsv "${ORIENTATION_TSV}" \
    --chr-count "${CHR_COUNT}"

require_file "${COPY2_ORIENTED}" "copy2 oriented FASTA"

remove_if_forced "${ALN_SAM}" "${ALN_BAM}" "${ALN_BAM}.bai" "${SYRI_OUT}" "${WORK_DIR}/${SYRI_PREFIX}syri.vcf" \
    "${WORK_DIR}/${SYRI_PREFIX}syri.summary" "${WORK_DIR}/pzi_copy1_vs_copy2.clean.png" "${WORK_DIR}/pzi_copy1_vs_copy2.detail.png"

if [ ! -s "${ALN_BAM}" ]; then
    echo ">>> Minimap2 SAM/BAM: oriented copy2 against copy1"
    minimap2 -ax "${MINIMAP_PRESET}" --eqx -t "${THREADS}" "${COPY1_FASTA}" "${COPY2_ORIENTED}" > "${ALN_SAM}"
    samtools sort -O BAM -@ "${THREADS}" "${ALN_SAM}" > "${ALN_BAM}"
    samtools index "${ALN_BAM}"
    rm -f "${ALN_SAM}"
else
    echo ">>> Existing BAM found, skipping: ${ALN_BAM}"
fi

if [ ! -s "${SYRI_OUT}" ]; then
    echo ">>> Running SyRI"
    syri -c "${ALN_BAM}" -r "${COPY1_FASTA}" -q "${COPY2_ORIENTED}" \
        -F B --nc "${THREADS}" --dir "${WORK_DIR}" --prefix "${SYRI_PREFIX}"
else
    echo ">>> Existing SyRI output found, skipping: ${SYRI_OUT}"
fi

require_file "${SYRI_OUT}" "copy1-vs-copy2 SyRI output"

{
    printf "%s\tcopy1\tlw:1.5\n" "${COPY1_FASTA}"
    printf "%s\tcopy2_oriented\tlw:1.5\n" "${COPY2_ORIENTED}"
} > "${GENOMES_TXT}"

if command -v plotsr >/dev/null 2>&1; then
    echo ">>> Running plotsr"
    plotsr --sr "${SYRI_OUT}" --genomes "${GENOMES_TXT}" \
        -o "${WORK_DIR}/pzi_copy1_vs_copy2.clean.png" -H 8 -W 6 -s 50000
    plotsr --sr "${SYRI_OUT}" --genomes "${GENOMES_TXT}" \
        -o "${WORK_DIR}/pzi_copy1_vs_copy2.detail.png" -H 8 -W 6
else
    echo "WARNING: plotsr not found in PATH; skipping plots." >&2
fi

{
    printf "copy1 vs copy2 haplotype synteny review\n"
    printf "=======================================\n\n"
    printf "Main checks:\n"
    printf -- "- Review %s: same_chr_aligned_bp should dominate over other_chr_aligned_bp.\n" "${SELF_SUMMARY}"
    printf -- "- Review %s: most chromosomes should have a clear plus/minus orientation assignment.\n" "${ORIENTATION_TSV}"
    printf -- "- Review %s and SyRI summary/VCF for large structural differences between haplotypes.\n" "${SYRI_OUT}"
    printf -- "- If plotsr outputs exist, use pzi_copy1_vs_copy2.clean.png for overview and detail.png for events.\n"
} > "${WORK_DIR}/review_checklist.txt"

echo
echo "=== copy1 vs copy2 SyRI outputs ==="
echo "Raw PAF                : ${SELF_PAF}"
echo "Raw alignment summary  : ${SELF_SUMMARY}"
echo "Orientation table      : ${ORIENTATION_TSV}"
echo "Oriented copy2 FASTA   : ${COPY2_ORIENTED}"
echo "BAM                    : ${ALN_BAM}"
echo "SyRI output            : ${SYRI_OUT}"
echo "Review checklist       : ${WORK_DIR}/review_checklist.txt"
