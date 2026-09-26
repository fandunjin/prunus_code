#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
COPY2_FASTA="${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy2.hap_provenance.chr16.fa"

WORK_DIR="${ROOT_DIR}/final_assembly/minimap_copy1_copy2"
THREADS=24
MINIMAP_PRESET="asm5"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
MINIMAP_ENV="minimap"
DOTPLOT_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py"
HEATMAP_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_heatmap.py"

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

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${COPY1_FASTA}" "P. zippeliana copy1 final n=16 FASTA"
require_file "${COPY2_FASTA}" "P. zippeliana copy2 final n=16 FASTA"
require_file "${DOTPLOT_SCRIPT}" "minimap synteny dotplot script"
require_file "${HEATMAP_SCRIPT}" "minimap synteny heatmap script"

mkdir -p "${WORK_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

prefix="${COPY1_NAME}_vs_${COPY2_NAME}_chr16"
paf="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.paf"
summary="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.ref_query_alignment_summary.tsv"
best="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.best_query_per_copy1_chr.tsv"
dotplot="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_dotplot.svg"
heatmap="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_heatmap.svg"

echo "=== PZI copy1 vs copy2 minimap2 synteny ==="
echo "copy1 reference: ${COPY1_FASTA}"
echo "copy2 query    : ${COPY2_FASTA}"
echo "Work dir       : ${WORK_DIR}"
echo "Preset         : ${MINIMAP_PRESET}"
echo "Threads        : ${THREADS}"
echo

remove_if_forced "${paf}" "${summary}" "${best}" "${dotplot}" "${heatmap}"

if [ ! -s "${paf}" ]; then
    echo ">>> Minimap2 PAF: ${COPY1_NAME} reference vs ${COPY2_NAME} query"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${COPY1_FASTA}" "${COPY2_FASTA}" > "${paf}"
else
    echo ">>> Existing PAF found, skipping: ${paf}"
fi

echo ">>> Summarizing copy1-copy2 PAF blocks"
summarize_paf "${paf}" "${summary}"
write_best_hits "${summary}" "${best}"

echo ">>> Drawing minimap2 synteny dotplot"
"${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
    --paf "${paf}" \
    --out-svg "${dotplot}" \
    --title "P. zippeliana copy1 vs copy2 minimap2 alignment"

echo ">>> Drawing minimap2 chromosome-pair heatmap"
"${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
    --summary "${summary}" \
    --out-svg "${heatmap}" \
    --title "P. zippeliana copy1 vs copy2 minimap2 aligned blocks" \
    --x-label "copy2 query chromosomes" \
    --y-label "copy1 reference chromosomes"

{
    printf "PZI copy1 vs copy2 minimap2 synteny review\n"
    printf "==========================================\n\n"
    printf "Main purpose:\n"
    printf -- "- Check whole-chromosome correspondence between the final copy1 and copy2 n=16 haplotype FASTAs.\n"
    printf -- "- This is a minimap2 PAF macro-synteny screen; use SyRI/plotsr outputs for detailed structural-variant interpretation.\n\n"
    printf "Main files:\n"
    printf -- "- %s\n" "${paf}"
    printf -- "- %s: copy1-by-copy2 total aligned block lengths.\n" "${summary}"
    printf -- "- %s: strongest copy2 chromosome hit for each copy1 chromosome.\n" "${best}"
    printf -- "- %s: coordinate-based minimap2 dotplot.\n" "${dotplot}"
    printf -- "- %s: chromosome-pair aligned-block heatmap.\n\n" "${heatmap}"
    printf "Interpretation notes:\n"
    printf -- "- Strong diagonal or one-to-one chromosome-pair signal supports the copy1/copy2 homolog pairing.\n"
    printf -- "- Large off-diagonal blocks should be checked against the final SyRI/plotsr and haplotype-provenance tables.\n"
} > "${WORK_DIR}/review_checklist.txt"

echo
echo "=== PZI copy1 vs copy2 minimap2 synteny outputs ==="
echo "Work dir        : ${WORK_DIR}"
echo "PAF             : ${paf}"
echo "Summary         : ${summary}"
echo "Best hits       : ${best}"
echo "Dotplot         : ${dotplot}"
echo "Heatmap         : ${heatmap}"
echo "Review checklist: ${WORK_DIR}/review_checklist.txt"
