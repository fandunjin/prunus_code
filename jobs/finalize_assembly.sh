#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
INPUT_FASTA="${ROOT_DIR}/haphic_n=16/merge/04.build/out_JBAT.FINAL.fa"

OUT_DIR="${ROOT_DIR}/final_assembly"
PREFIX="pzi_2n32"

# The first 32 records in out_JBAT.FINAL.fa are treated as chromosomes.
CHR_COUNT=32

# Homolog pairing used by HapHiC merge result after manual Juicebox curation:
# pair01 = chr01 + chr02, pair02 = chr03 + chr04, ..., pair16 = chr31 + chr32.
PAIR_MODE="adjacent"

# Representative n=16 copy to emit for reference comparison.
# 1 selects chr01,chr03,...,chr31; 2 selects chr02,chr04,...,chr32.
REPRESENTATIVE_COPY=1

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SAMTOOLS_ENV="haphic"
FINALIZE_SCRIPT="${ROOT_DIR}/scripts/finalize_assembly.py"

# ===================== Pipeline starts here =====================

FINAL_FASTA="${OUT_DIR}/${PREFIX}.final.fa"
CHR_FASTA="${OUT_DIR}/${PREFIX}.chromosomes.fa"
UNPLACED_FASTA="${OUT_DIR}/${PREFIX}.unplaced.fa"
REP_FASTA="${OUT_DIR}/${PREFIX}.n16_representative_copy${REPRESENTATIVE_COPY}.fa"

NAME_MAP="${OUT_DIR}/${PREFIX}.rename_map.tsv"
HOMOLOG_PAIRS="${OUT_DIR}/${PREFIX}.homolog_pairs.tsv"
CHROM_SIZES="${OUT_DIR}/${PREFIX}.chrom.sizes"
CHR_SIZES="${OUT_DIR}/${PREFIX}.chromosomes.sizes"
SEQ_STATS="${OUT_DIR}/${PREFIX}.sequence_stats.tsv"
GAP_STATS="${OUT_DIR}/${PREFIX}.gap_stats.tsv"
AGP_FILE="${OUT_DIR}/${PREFIX}.agp"
MD5_FILE="${OUT_DIR}/${PREFIX}.md5"

require_file() {
    local file_path="$1"
    local label="$2"

    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${INPUT_FASTA}" "Juicebox final FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${FINALIZE_SCRIPT}" "Finalize Python script"

mkdir -p "${OUT_DIR}"

source "${CONDA_SH}"
conda activate "${SAMTOOLS_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo "=== Finalizing 2n=32 chromosome assembly ==="
echo "Input FASTA : ${INPUT_FASTA}"
echo "Output dir  : ${OUT_DIR}"
echo "Prefix      : ${PREFIX}"
echo "Chr count   : ${CHR_COUNT}"
echo "Pair mode   : ${PAIR_MODE}"
echo

"${PYTHON_BIN}" "${FINALIZE_SCRIPT}" \
    --input-fasta "${INPUT_FASTA}" \
    --out-dir "${OUT_DIR}" \
    --prefix "${PREFIX}" \
    --chr-count "${CHR_COUNT}" \
    --pair-mode "${PAIR_MODE}" \
    --representative-copy "${REPRESENTATIVE_COPY}"

samtools faidx "${FINAL_FASTA}"
samtools faidx "${CHR_FASTA}"
samtools faidx "${REP_FASTA}"

echo
echo "=== Final assembly files ==="
echo "Final FASTA          : ${FINAL_FASTA}"
echo "Chromosome FASTA     : ${CHR_FASTA}"
echo "Unplaced FASTA       : ${UNPLACED_FASTA}"
echo "Representative FASTA : ${REP_FASTA}"
echo "Rename map           : ${NAME_MAP}"
echo "Homolog pairs        : ${HOMOLOG_PAIRS}"
echo "Chrom sizes          : ${CHROM_SIZES}"
echo "Sequence stats       : ${SEQ_STATS}"
echo "Gap stats            : ${GAP_STATS}"
echo "AGP                  : ${AGP_FILE}"
echo "MD5                  : ${MD5_FILE}"
