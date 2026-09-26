#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
HOMOLOG_PAIRS="${ROOT_DIR}/final_assembly/pzi_2n32.homolog_pairs.tsv"

# Original hifiasm haplotig contigs after FCS cleanup.
HAP1_FASTA="${ROOT_DIR}/fcs/final_clean/ONT.asm.HD-14-2DNA.hic.hap1.p_ctg.final_clean.fasta"
HAP2_FASTA="${ROOT_DIR}/fcs/final_clean/ONT.asm.HD-14-2DNA.hic.hap2.p_ctg.final_clean.fasta"

OUT_DIR="${ROOT_DIR}/final_assembly/haplotype_provenance"
PREFIX="pzi_2n32"

THREADS=24
MINIMAP_PRESET="asm5"
MIN_MAPQ=10
MIN_ALN_LEN=10000
DOMINANCE_THRESHOLD=0.60

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
MINIMAP_ENV="minimap"
PROVENANCE_SCRIPT="${ROOT_DIR}/scripts/haplotype_provenance.py"

# ===================== Pipeline starts here =====================

HAP1_PAF="${OUT_DIR}/${PREFIX}.hap1_to_final.paf"
HAP2_PAF="${OUT_DIR}/${PREFIX}.hap2_to_final.paf"
OUT_PREFIX="${OUT_DIR}/${PREFIX}"
COPY1_FASTA="${OUT_DIR}/${PREFIX}.copy1.hap_provenance.chr16.fa"
COPY2_FASTA="${OUT_DIR}/${PREFIX}.copy2.hap_provenance.chr16.fa"

require_file() {
    local file_path="$1"
    local label="$2"

    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${FINAL_FASTA}" "Final FASTA"
require_file "${HOMOLOG_PAIRS}" "Homolog pairs table"
require_file "${HAP1_FASTA}" "hifiasm hap1 FASTA"
require_file "${HAP2_FASTA}" "hifiasm hap2 FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${PROVENANCE_SCRIPT}" "haplotype provenance Python script"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo "=== Final haplotype provenance check ==="
echo "Final FASTA : ${FINAL_FASTA}"
echo "Hap1 FASTA  : ${HAP1_FASTA}"
echo "Hap2 FASTA  : ${HAP2_FASTA}"
echo "Output dir  : ${OUT_DIR}"
echo

if [ ! -s "${HAP1_PAF}" ]; then
    echo ">>> Aligning hifiasm hap1 contigs to final assembly"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${FINAL_FASTA}" "${HAP1_FASTA}" > "${HAP1_PAF}"
else
    echo ">>> ${HAP1_PAF} exists, skipping hap1 alignment"
fi

if [ ! -s "${HAP2_PAF}" ]; then
    echo ">>> Aligning hifiasm hap2 contigs to final assembly"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${FINAL_FASTA}" "${HAP2_FASTA}" > "${HAP2_PAF}"
else
    echo ">>> ${HAP2_PAF} exists, skipping hap2 alignment"
fi

echo ">>> Inferring hap1/hap2 provenance and writing corrected copy FASTAs"
"${PYTHON_BIN}" "${PROVENANCE_SCRIPT}" \
    --final-fasta "${FINAL_FASTA}" \
    --hap1-paf "${HAP1_PAF}" \
    --hap2-paf "${HAP2_PAF}" \
    --homolog-pairs "${HOMOLOG_PAIRS}" \
    --out-prefix "${OUT_PREFIX}" \
    --copy1-fasta "${COPY1_FASTA}" \
    --copy2-fasta "${COPY2_FASTA}" \
    --min-mapq "${MIN_MAPQ}" \
    --min-aln-len "${MIN_ALN_LEN}" \
    --dominance-threshold "${DOMINANCE_THRESHOLD}"

samtools faidx "${COPY1_FASTA}"
samtools faidx "${COPY2_FASTA}"

{
    printf "Final haplotype provenance review\n"
    printf "=================================\n\n"
    printf "Required checks:\n"
    printf -- "- Review %s.chr_hap_provenance.tsv for chromosomes marked ambiguous.\n" "${OUT_PREFIX}"
    printf -- "- Review %s.corrected_homolog_pairs.tsv for swapped or ambiguous_keep_original pairs.\n" "${OUT_PREFIX}"
    printf -- "- Use the corrected copy FASTAs for downstream SyRI/plotsr when all important chromosomes are confident.\n"
    printf -- "- Pay special attention to final chr27-chr30 / normalized pair14-pair15.\n"
} > "${OUT_DIR}/review_checklist.txt"

echo
echo "=== Haplotype provenance outputs ==="
echo "Hap1 PAF             : ${HAP1_PAF}"
echo "Hap2 PAF             : ${HAP2_PAF}"
echo "Chromosome summary   : ${OUT_PREFIX}.chr_hap_provenance.tsv"
echo "Corrected pairs      : ${OUT_PREFIX}.corrected_homolog_pairs.tsv"
echo "Corrected copy1 FASTA: ${COPY1_FASTA}"
echo "Corrected copy2 FASTA: ${COPY2_FASTA}"
