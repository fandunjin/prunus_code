#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

OUT_DIR="${ROOT_DIR}/final_assembly/qc"
THREADS=48

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
QUAST_ENV="quast"
BUSCO_ENV="busco"
FCS_ENV="fcs"
MAP_ENV="minimap"
FIND_FIRST_SCRIPT="${ROOT_DIR}/scripts/find_first_file.py"

BUSCO_LINEAGE="${HOME}/data/busco_db/eudicotyledons_odb12/"

# FCS-GX settings.
RUN_FCS_GX=1
TAX_ID="3754"
GX_DB="${HOME}/cx/docker_images/fcs/fcs_db"
GX_IMG="${HOME}/cx/docker_images/fcs/fcs-gx.sif"
ADAPTOR_IMG="${HOME}/cx/docker_images/fcs/fcs-adaptor.sif"
BIND_PATHS="-B ${PROJECT_ROOT}:${PROJECT_ROOT}"

# Optional long-read depth check. Set ONT_READS to the final ONT FASTQ if available.
RUN_ONT_DEPTH=0
ONT_READS=""

# Optional Merqury check. Set both variables if meryl db and merqury are available.
RUN_MERQURY=0
MERYL_DB=""

# ===================== Pipeline starts here =====================

require_file() {
    local file_path="$1"
    local label="$2"

    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_path() {
    local file_path="$1"
    local label="$2"

    if [ ! -e "${file_path}" ]; then
        echo "ERROR: ${label} does not exist: ${file_path}" >&2
        exit 1
    fi
}

require_file "${FINAL_FASTA}" "Final FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${FIND_FIRST_SCRIPT}" "find-first-file Python script"

mkdir -p "${OUT_DIR}/quast" "${OUT_DIR}/busco" "${OUT_DIR}/fcs" "${OUT_DIR}/read_depth" "${OUT_DIR}/merqury"

echo "=== Final assembly QC ==="
echo "Final FASTA : ${FINAL_FASTA}"
echo "Output dir  : ${OUT_DIR}"
echo

source "${CONDA_SH}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo ">>> QUAST"
conda activate "${QUAST_ENV}"
quast -t "${THREADS}" -o "${OUT_DIR}/quast" "${FINAL_FASTA}"

echo ">>> BUSCO"
conda activate "${BUSCO_ENV}"
busco -i "${FINAL_FASTA}" -m genome \
    -l "${BUSCO_LINEAGE}" \
    -c "${THREADS}" --offline \
    -o "busco_pzi_2n32_final" \
    --out_path "${OUT_DIR}/busco" \
    --force

if [ "${RUN_FCS_GX}" -eq 1 ]; then
    echo ">>> FCS adaptor + GX contamination screen"
    conda activate "${FCS_ENV}"
    mkdir -p "${OUT_DIR}/fcs/adaptor_out" "${OUT_DIR}/fcs/gx_out" "${OUT_DIR}/fcs/final_clean"

    apptainer exec ${BIND_PATHS} "${ADAPTOR_IMG}" \
        /app/fcs/bin/av_screen_x \
        -o "${OUT_DIR}/fcs/adaptor_out/" \
        --euk "${FINAL_FASTA}"

    ADAPTOR_CLEAN=$("${PYTHON_BIN}" "${FIND_FIRST_SCRIPT}" \
        --dir "${OUT_DIR}/fcs/adaptor_out/cleaned_sequences" \
        --pattern "*")
    require_file "${ADAPTOR_CLEAN}" "FCS adaptor-cleaned FASTA"

    apptainer exec ${BIND_PATHS} "${GX_IMG}" \
        python3 /app/bin/run_gx \
        --fasta "${ADAPTOR_CLEAN}" \
        --out-dir "${OUT_DIR}/fcs/gx_out/" \
        --gx-db "${GX_DB}" \
        --tax-id "${TAX_ID}"

    GX_REPORT=$("${PYTHON_BIN}" "${FIND_FIRST_SCRIPT}" \
        --dir "${OUT_DIR}/fcs/gx_out" \
        --pattern "*fcs_gx_report.txt")
    require_file "${GX_REPORT}" "FCS-GX report"

    apptainer exec ${BIND_PATHS} "${GX_IMG}" \
        /app/bin/gx clean-genome \
        --input "${ADAPTOR_CLEAN}" \
        --action-report "${GX_REPORT}" \
        --output "${OUT_DIR}/fcs/final_clean/pzi_2n32.final.fcs_clean.fasta" \
        --contam-fasta-out "${OUT_DIR}/fcs/final_clean/pzi_2n32.final.contam.fasta"
else
    echo ">>> RUN_FCS_GX=0, skipping FCS contamination screen"
fi

if [ "${RUN_ONT_DEPTH}" -eq 1 ]; then
    require_file "${ONT_READS}" "ONT reads"
    echo ">>> ONT read depth check"
    conda activate "${MAP_ENV}"
    minimap2 -ax map-ont -t "${THREADS}" "${FINAL_FASTA}" "${ONT_READS}" \
        | samtools sort -@ "${THREADS}" -o "${OUT_DIR}/read_depth/ont_to_final.sorted.bam"
    samtools index -@ "${THREADS}" "${OUT_DIR}/read_depth/ont_to_final.sorted.bam"
    samtools depth -aa "${OUT_DIR}/read_depth/ont_to_final.sorted.bam" \
        | awk 'BEGIN{OFS="\t"} {sum[$1]+=$3; n[$1]++} END{print "seq_id","mean_depth"; for (s in sum) print s,sum[s]/n[s]}' \
        > "${OUT_DIR}/read_depth/ont_mean_depth_by_seq.tsv"
else
    echo ">>> RUN_ONT_DEPTH=0, skipping ONT depth check"
fi

if [ "${RUN_MERQURY}" -eq 1 ]; then
    require_path "${MERYL_DB}" "Merqury/Meryl database"
    echo ">>> Merqury"
    merqury.sh "${MERYL_DB}" "${FINAL_FASTA}" "${OUT_DIR}/merqury/pzi_2n32_final"
else
    echo ">>> RUN_MERQURY=0, skipping Merqury"
fi

{
    printf "Final assembly QC review\n"
    printf "========================\n\n"
    printf "Required checks:\n"
    printf -- "- Compare QUAST metrics against previous hifiasm/HapHiC versions.\n"
    printf -- "- Compare BUSCO complete/duplicated/missing values against previous versions.\n"
    printf -- "- Review FCS-GX report and contaminant FASTA, if RUN_FCS_GX=1.\n"
    printf -- "- If ONT depth or Merqury was enabled, flag chromosomes with abnormal depth/QV/completeness.\n"
} > "${OUT_DIR}/qc_summary.todo.txt"

echo
echo "=== Final assembly QC outputs ==="
echo "QUAST      : ${OUT_DIR}/quast"
echo "BUSCO      : ${OUT_DIR}/busco"
echo "FCS        : ${OUT_DIR}/fcs"
echo "Read depth : ${OUT_DIR}/read_depth"
echo "Merqury    : ${OUT_DIR}/merqury"
