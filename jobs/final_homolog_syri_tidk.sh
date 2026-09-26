#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
HOMOLOG_PAIRS="${ROOT_DIR}/final_assembly/pzi_2n32.homolog_pairs.tsv"
CORRECTED_HOMOLOG_PAIRS="${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv"

REF_NAME="ppa"
REF_LABEL="ppa"
REF_GENOME="${ROOT_DIR}/data/ppa/GCA_964059325.1_drPruPadu1.1_alternate_haplotype_genomic.fna"
REF_CHR_MODE="ppa_super"

QUERY_NAME="pzi"
QUERY_LABEL="pzi"

WORK_DIR="${ROOT_DIR}/final_assembly/homolog_syri_tidk"

THREADS=24
MINIMAP_PRESET="asm5"
SYRI_CHR_COUNT=16

TELO_REPEAT="TTTAGGG"
TIDK_WINDOW=5000

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SYRI_ENV="syri_run"
TIDK_ENV="haphic"
HOMOLOG_HELPER_SCRIPT="${ROOT_DIR}/scripts/homolog_syri_helpers.py"

FORCE_RERUN=0

# ===================== Pipeline starts here =====================

REF_RENAMED="${WORK_DIR}/${REF_NAME}.chr${SYRI_CHR_COUNT}.renamed.fa"
COPY1_FASTA="${WORK_DIR}/${QUERY_NAME}.copy1.chr${SYRI_CHR_COUNT}.renamed.fa"
COPY2_FASTA="${WORK_DIR}/${QUERY_NAME}.copy2.chr${SYRI_CHR_COUNT}.renamed.fa"
COPY_MAP="${WORK_DIR}/${QUERY_NAME}.copy_fastas.map.tsv"
SELF_PAF="${WORK_DIR}/${QUERY_NAME}.copy1_vs_copy2.paf"
SELF_SUMMARY="${WORK_DIR}/${QUERY_NAME}.copy1_vs_copy2.same_chr_alignment.tsv"

TIDK_PREFIX="tidk_${QUERY_NAME}_2n32_final"

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

require_file "${FINAL_FASTA}" "Final FASTA"
require_file "${REF_GENOME}" "Reference genome"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${HOMOLOG_HELPER_SCRIPT}" "Homolog/SyRI helper Python script"

if [ -s "${CORRECTED_HOMOLOG_PAIRS}" ]; then
    ACTIVE_HOMOLOG_PAIRS="${CORRECTED_HOMOLOG_PAIRS}"
else
    ACTIVE_HOMOLOG_PAIRS="${HOMOLOG_PAIRS}"
fi
require_file "${ACTIVE_HOMOLOG_PAIRS}" "Active homolog pairs table"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

source "${CONDA_SH}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo "=== Final homolog, telomere, and SyRI/TIDK analysis ==="
echo "Final FASTA : ${FINAL_FASTA}"
echo "Pairs table : ${ACTIVE_HOMOLOG_PAIRS}"
echo "Reference   : ${REF_GENOME}"
echo "Work dir    : ${WORK_DIR}"
echo

echo ">>> Building copy1/copy2 n=16 FASTAs and normalized reference"
"${PYTHON_BIN}" "${HOMOLOG_HELPER_SCRIPT}" build-fastas \
    --final-fasta "${FINAL_FASTA}" \
    --homolog-pairs "${ACTIVE_HOMOLOG_PAIRS}" \
    --ref-fasta "${REF_GENOME}" \
    --ref-out "${REF_RENAMED}" \
    --copy1-out "${COPY1_FASTA}" \
    --copy2-out "${COPY2_FASTA}" \
    --copy-map "${COPY_MAP}" \
    --ref-chr-mode "${REF_CHR_MODE}" \
    --chr-count "${SYRI_CHR_COUNT}"

conda activate "${SYRI_ENV}"

remove_if_forced "${SELF_PAF}" "${SELF_SUMMARY}"
if [ ! -s "${SELF_PAF}" ]; then
    echo ">>> Aligning copy1 against copy2 for homolog-pair support"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${COPY1_FASTA}" "${COPY2_FASTA}" > "${SELF_PAF}"
else
    echo ">>> ${SELF_PAF} exists, skipping copy1-vs-copy2 alignment"
fi

"${PYTHON_BIN}" "${HOMOLOG_HELPER_SCRIPT}" summarize-copy-alignment \
    --paf "${SELF_PAF}" \
    --out-tsv "${SELF_SUMMARY}" \
    --chr-count "${SYRI_CHR_COUNT}"

echo ">>> TIDK telomere search on final 2n=32 FASTA"
conda activate "${TIDK_ENV}"
remove_if_forced "${TIDK_PREFIX}_telomeric_repeat_windows.tsv" "${TIDK_PREFIX}_plot.svg"
if [ ! -s "${TIDK_PREFIX}_telomeric_repeat_windows.tsv" ]; then
    tidk search \
        --string "${TELO_REPEAT}" \
        --window "${TIDK_WINDOW}" \
        --output "${TIDK_PREFIX}" \
        --dir "${WORK_DIR}" \
        "${FINAL_FASTA}"
else
    echo ">>> ${TIDK_PREFIX}_telomeric_repeat_windows.tsv exists, skipping tidk search"
fi
tidk plot \
    --tsv "${TIDK_PREFIX}_telomeric_repeat_windows.tsv" \
    --output "${TIDK_PREFIX}_plot"

orient_query_fasta() {
    local query_in="$1"
    local query_out="$2"
    local prefix="$3"
    local orient_paf="${prefix}.orientation.paf"
    local orient_tsv="${prefix}.orientation.tsv"

    if [ -s "${query_out}" ] && [ "${FORCE_RERUN}" -ne 1 ]; then
        echo "    -> ${query_out} exists, skipping orientation"
        return
    fi

    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${REF_RENAMED}" "${query_in}" > "${orient_paf}"

    "${PYTHON_BIN}" "${HOMOLOG_HELPER_SCRIPT}" orient-query \
        --paf "${orient_paf}" \
        --orientation-tsv "${orient_tsv}" \
        --query-in "${query_in}" \
        --query-out "${query_out}" \
        --chr-count "${SYRI_CHR_COUNT}"
}

run_syri_copy() {
    local copy_id="$1"
    local query_in="$2"
    local prefix="${REF_NAME}_vs_${QUERY_NAME}_${copy_id}_chr${SYRI_CHR_COUNT}"
    local query_oriented="${WORK_DIR}/${QUERY_NAME}.${copy_id}.chr${SYRI_CHR_COUNT}.oriented.fa"
    local aln_sam="${prefix}.sam"
    local aln_bam="${prefix}.bam"
    local syri_out="${prefix}.syri.out"
    local genomes_txt="${prefix}.genomes.txt"

    echo ">>> SyRI/plotsr for ${copy_id}"
    conda activate "${SYRI_ENV}"
    orient_query_fasta "${query_in}" "${query_oriented}" "${prefix}"

    remove_if_forced "${aln_sam}" "${aln_bam}" "${aln_bam}.bai" "${syri_out}" "${genomes_txt}" \
        "${prefix}_clean.png" "${prefix}_detail.png"

    if [ ! -s "${aln_bam}" ]; then
        minimap2 -ax "${MINIMAP_PRESET}" --eqx -t "${THREADS}" "${REF_RENAMED}" "${query_oriented}" > "${aln_sam}"
        samtools sort -O BAM -@ "${THREADS}" "${aln_sam}" > "${aln_bam}"
        samtools index "${aln_bam}"
        rm -f "${aln_sam}"
    else
        echo "    -> ${aln_bam} exists, skipping minimap2"
    fi

    if [ ! -s "${syri_out}" ]; then
        syri -c "${aln_bam}" -r "${REF_RENAMED}" -q "${query_oriented}" \
            -F B --nc "${THREADS}" --prefix "${prefix}."
    else
        echo "    -> ${syri_out} exists, skipping SyRI"
    fi
    require_file "${syri_out}" "SyRI output for ${copy_id}"

    {
        printf "%s\t%s\tlw:1.5\n" "${REF_RENAMED}" "${REF_LABEL}"
        printf "%s\t%s_%s\tlw:1.5\n" "${query_oriented}" "${QUERY_LABEL}" "${copy_id}"
    } > "${genomes_txt}"

    plotsr --sr "${syri_out}" --genomes "${genomes_txt}" \
        -o "${prefix}_clean.png" -H 8 -W 6 -s 50000
    plotsr --sr "${syri_out}" --genomes "${genomes_txt}" \
        -o "${prefix}_detail.png" -H 8 -W 6
}

run_syri_copy "copy1" "${COPY1_FASTA}"
run_syri_copy "copy2" "${COPY2_FASTA}"

{
    printf "Final homolog/SyRI/TIDK review\n"
    printf "==============================\n\n"
    printf "Required checks:\n"
    printf -- "- Review %s; each homolog pair should align mostly to the same normalized chromosome.\n" "${SELF_SUMMARY}"
    printf -- "- Review %s_plot.svg; telomeric repeat peaks should be concentrated at chromosome ends.\n" "${TIDK_PREFIX}"
    printf -- "- Review ppa_vs_pzi_copy1_chr16_clean.png and ppa_vs_pzi_copy2_chr16_clean.png for large-scale collinearity.\n"
    printf -- "- Record any real inversions/translocations separately from likely assembly errors.\n"
} > "${WORK_DIR}/review_checklist.txt"

echo
echo "=== Final homolog/SyRI/TIDK outputs ==="
echo "Copy FASTA map        : ${COPY_MAP}"
echo "Copy1 FASTA           : ${COPY1_FASTA}"
echo "Copy2 FASTA           : ${COPY2_FASTA}"
echo "Copy alignment PAF    : ${SELF_PAF}"
echo "Copy alignment summary: ${SELF_SUMMARY}"
echo "TIDK plot             : ${WORK_DIR}/${TIDK_PREFIX}_plot.svg"
echo "SyRI/plotsr outputs   : ${WORK_DIR}"
