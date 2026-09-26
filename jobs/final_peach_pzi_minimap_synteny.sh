#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Final Prunus zippeliana haplotype-provenance n=16 assemblies.
# Use these for x=8 peach vs n=16 P. zippeliana macro-synteny checks.
PZI_NAME="pzi"
PZI_COPY1_FASTA="${ROOT_DIR}/final_genome/pzi_copy1.fa"
PZI_COPY2_FASTA="${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy2.hap_provenance.chr16.fa"

WORK_DIR="${ROOT_DIR}/final_assembly/minimap_pzi_peach"
THREADS=16
MINIMAP_PRESET="asm5"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
MINIMAP_ENV="minimap"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py"

FORCE_RERUN="${FORCE_RERUN:-0}"

# ===================== Pipeline starts here =====================

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

summarize_paf() {
    local paf="$1"
    local out_tsv="$2"
    local tmp_tsv="${out_tsv}.tmp.$$"

    require_file "${paf}" "PAF alignment"

    awk 'BEGIN {
             OFS="\t";
         }
         {
             key=$6 SUBSEP $1;
             aligned[key]+=$11;
             matching[key]+=$10;
             records[key]+=1;
         }
         END {
             for (key in aligned) {
                 split(key, a, SUBSEP);
                 print a[1], a[2], aligned[key], matching[key], records[key];
             }
         }' "${paf}" | sort -k1,1 -k3,3nr > "${tmp_tsv}"

    {
        printf "reference_chr\tquery_chr\taligned_block_bp\tmatching_bp\talignment_records\n"
        cat "${tmp_tsv}"
    } > "${out_tsv}"
    rm -f "${tmp_tsv}"
}

write_best_hits() {
    local summary_tsv="$1"
    local out_tsv="$2"

    require_file "${summary_tsv}" "PAF summary table"

    awk 'BEGIN {
             OFS="\t";
             print "reference_chr","best_query_chr","aligned_block_bp","matching_bp","alignment_records";
         }
         NR == 1 { next }
         !seen[$1] {
             print $1, $2, $3, $4, $5;
             seen[$1]=1;
         }' "${summary_tsv}" > "${out_tsv}"
}

run_copy_alignment() {
    local copy_id="$1"
    local query_fasta="$2"

    local prefix="${PEACH_NAME}_vs_${PZI_NAME}_${copy_id}_chr16"
    local paf="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.paf"
    local summary="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.ref_query_alignment_summary.tsv"
    local best="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.best_query_per_peach_chr.tsv"
    local dotplot="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_dotplot.svg"

    require_file "${query_fasta}" "P. zippeliana ${copy_id} final n=16 FASTA"
    require_file "${PLOT_SCRIPT}" "minimap synteny dotplot script"

    remove_if_forced "${paf}" "${summary}" "${best}" "${dotplot}"

    if [ ! -s "${paf}" ]; then
        echo ">>> Minimap2 PAF: ${PEACH_NAME} reference vs ${PZI_NAME} ${copy_id}"
        minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PEACH_FASTA}" "${query_fasta}" > "${paf}"
    else
        echo ">>> Existing PAF found, skipping: ${paf}"
    fi

    echo ">>> Summarizing PAF blocks for ${copy_id}"
    summarize_paf "${paf}" "${summary}"
    write_best_hits "${summary}" "${best}"

    echo ">>> Drawing minimap2 synteny dotplot for ${copy_id}"
    "${PYTHON_BIN}" "${PLOT_SCRIPT}" \
        --paf "${paf}" \
        --out-svg "${dotplot}" \
        --title "Peach vs P. zippeliana ${copy_id} minimap2 alignment"

    echo "    PAF        : ${paf}"
    echo "    Summary    : ${summary}"
    echo "    Best hits  : ${best}"
    echo "    Dotplot    : ${dotplot}"
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${PEACH_FASTA}" "Peach reference FASTA"
require_file "${PLOT_SCRIPT}" "minimap synteny dotplot script"

mkdir -p "${WORK_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Peach vs Prunus zippeliana final assembly minimap2 synteny ==="
echo "Peach reference : ${PEACH_FASTA}"
echo "PZI copy1 FASTA : ${PZI_COPY1_FASTA}"
echo "PZI copy2 FASTA : ${PZI_COPY2_FASTA}"
echo "Work dir        : ${WORK_DIR}"
echo "Preset          : ${MINIMAP_PRESET}"
echo "Threads         : ${THREADS}"
echo

run_copy_alignment "copy1" "${PZI_COPY1_FASTA}"
run_copy_alignment "copy2" "${PZI_COPY2_FASTA}"

{
    printf "Peach vs P. zippeliana minimap2 synteny review\n"
    printf "==============================================\n\n"
    printf "Main purpose:\n"
    printf -- "- Test whether the 8 peach chromosomes show an approximate 1-to-2 macro-syntenic relationship with the 16 P. zippeliana haploid chromosomes.\n"
    printf -- "- This is a genome-scale alignment screen. It does not by itself prove WGD; gene-level collinearity and Ks analyses are still required.\n\n"
    printf "Main files:\n"
    printf -- "- %s/%s_vs_%s_copy1_chr16.%s.paf\n" "${WORK_DIR}" "${PEACH_NAME}" "${PZI_NAME}" "${MINIMAP_PRESET}"
    printf -- "- %s/%s_vs_%s_copy2_chr16.%s.paf\n" "${WORK_DIR}" "${PEACH_NAME}" "${PZI_NAME}" "${MINIMAP_PRESET}"
    printf -- "- *.ref_query_alignment_summary.tsv: all peach-chromosome by PZI-chromosome alignment totals.\n"
    printf -- "- *.best_query_per_peach_chr.tsv: strongest PZI chromosome hit for each peach chromosome.\n\n"
    printf -- "- *.alignment_dotplot.svg: coordinate-based minimap2 synteny dotplots.\n\n"
    printf "Interpretation notes:\n"
    printf -- "- A repeated pattern in which one peach chromosome aligns to two PZI chromosomes is compatible with WGD, but also with chromosome fission/dysploidy.\n"
    printf -- "- Follow up with MCScan/JCVI syntenic-depth and Ks analyses before writing WGD as a conclusion.\n"
} > "${WORK_DIR}/review_checklist.txt"

echo
echo "=== Peach vs PZI minimap2 synteny outputs ==="
echo "Work dir        : ${WORK_DIR}"
echo "Review checklist: ${WORK_DIR}/review_checklist.txt"
