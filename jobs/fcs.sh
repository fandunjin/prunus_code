#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate fcs

work_path="${PROJECT_ROOT}/fcs"
tax_id="3754"

input_dir="${PROJECT_ROOT}/hifiasm_output"

name1="ONT.asm.HD-14-2DNA.hic.hap1.p_ctg.fasta"
name2="ONT.asm.HD-14-2DNA.hic.hap2.p_ctg.fasta"
p_ctg="ONT.asm.HD-14-2DNA.hic.p_ctg.fasta"

fasta_hap1="${input_dir}/${name1}"
fasta_hap2="${input_dir}/${name2}"
fasta_p_ctg="${input_dir}/${p_ctg}"

gx_db="${HOME}/cx/docker_images/fcs/fcs_db"
gx_img="${HOME}/cx/docker_images/fcs/fcs-gx.sif"
adaptor_img="${HOME}/cx/docker_images/fcs/fcs-adaptor.sif"

BIND_PATHS="-B ${PROJECT_ROOT}:${PROJECT_ROOT}"

cd ${work_path}

mkdir -p ${work_path}/adaptor_out_hap1 ${work_path}/adaptor_out_hap2 ${work_path}/adaptor_out_p_ctg
mkdir -p ${work_path}/gx_out_hap1 ${work_path}/gx_out_hap2 ${work_path}/gx_out_p_ctg
mkdir -p ${work_path}/final_clean

echo "Step 1: Running FCS-Adaptor..."
apptainer exec ${BIND_PATHS} "${adaptor_img}" /app/fcs/bin/av_screen_x -o ${work_path}/adaptor_out_hap1/ --euk "${fasta_hap1}"
apptainer exec ${BIND_PATHS} "${adaptor_img}" /app/fcs/bin/av_screen_x -o ${work_path}/adaptor_out_hap2/ --euk "${fasta_hap2}"
apptainer exec ${BIND_PATHS} "${adaptor_img}" /app/fcs/bin/av_screen_x -o ${work_path}/adaptor_out_p_ctg/ --euk "${fasta_p_ctg}"

adapt_clean_hap1="${work_path}/adaptor_out_hap1/cleaned_sequences/${name1}"
adapt_clean_hap2="${work_path}/adaptor_out_hap2/cleaned_sequences/${name2}"
adapt_clean_p_ctg="${work_path}/adaptor_out_p_ctg/cleaned_sequences/${p_ctg}"

echo "Step 2: Running FCS-GX (This will take hours and high memory)..."
apptainer exec ${BIND_PATHS} "${gx_img}" python3 /app/bin/run_gx --fasta "${adapt_clean_hap1}" --out-dir ${work_path}/gx_out_hap1/ --gx-db "${gx_db}" --tax-id ${tax_id}
apptainer exec ${BIND_PATHS} "${gx_img}" python3 /app/bin/run_gx --fasta "${adapt_clean_hap2}" --out-dir ${work_path}/gx_out_hap2/ --gx-db "${gx_db}" --tax-id ${tax_id}
apptainer exec ${BIND_PATHS} "${gx_img}" python3 /app/bin/run_gx --fasta "${adapt_clean_p_ctg}" --out-dir ${work_path}/gx_out_p_ctg/ --gx-db "${gx_db}" --tax-id ${tax_id}

run_clean() {
    local INPUT_FASTA=$1
    local GX_DIR=$2
    local FINAL_DIR=$3
    local NAME=$(basename "$INPUT_FASTA")
    
    local GX_REPORT=$(find "${GX_DIR}" -name "*fcs_gx_report.txt" | head -n 1)

    if [[ -z "$GX_REPORT" ]]; then
        echo "Error: GX Report not found in ${GX_DIR}. Check fcs_err log for OOM or path errors."
        return
    fi

    echo "Cleaning ${NAME}..."
    apptainer exec ${BIND_PATHS} "${gx_img}" /app/bin/gx clean-genome \
        --input "$INPUT_FASTA" \
        --action-report "$GX_REPORT" \
        --output "${FINAL_DIR}/${NAME%.*}.final_clean.fasta" \
        --contam-fasta-out "${FINAL_DIR}/${NAME}.contam.fasta"
}

echo "Step 3: Starting final cleanup..."
run_clean "${adapt_clean_hap1}" "${work_path}/gx_out_hap1" "${work_path}/final_clean"
run_clean "${adapt_clean_hap2}" "${work_path}/gx_out_hap2" "${work_path}/final_clean"
run_clean "${adapt_clean_p_ctg}" "${work_path}/gx_out_p_ctg" "${work_path}/final_clean"

echo "FCS Pipeline Finished."