#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FINAL_DIR="${ROOT_DIR}/final"
PREFIX="pzi_2n32"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SAMTOOLS_ENV="haphic"
PACKAGE_SCRIPT="${ROOT_DIR}/scripts/package_final_assembly.py"

# Set to 1 only if you want to copy very large raw evidence files
# such as BAM/VCF/SyRI out files into final/qc/syri_plotsr/large_evidence.
INCLUDE_LARGE_EVIDENCE=0

# ===================== Pipeline starts here =====================

FINAL_FASTA="${SOURCE_DIR}/${PREFIX}.final.fa"
CHR_FASTA="${SOURCE_DIR}/${PREFIX}.chromosomes.fa"
UNPLACED_FASTA="${SOURCE_DIR}/${PREFIX}.unplaced.fa"
COPY1_FASTA="${SOURCE_DIR}/haplotype_provenance/${PREFIX}.copy1.hap_provenance.chr16.fa"
COPY2_FASTA="${SOURCE_DIR}/haplotype_provenance/${PREFIX}.copy2.hap_provenance.chr16.fa"
COPY1_FASTA_ALT="${SOURCE_DIR}/${PREFIX}.copy1.hap_provenance.chr16.fa"
COPY2_FASTA_ALT="${SOURCE_DIR}/${PREFIX}.copy2.hap_provenance.chr16.fa"

require_file() {
    local file_path="$1"
    local label="$2"

    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${FINAL_FASTA}" "Final 2n=32 FASTA"
require_file "${CHR_FASTA}" "Chromosome-only FASTA"
require_file "${UNPLACED_FASTA}" "Unplaced FASTA"
require_file "${PACKAGE_SCRIPT}" "Final package Python script"
require_file "${CONDA_SH}" "Conda initialization script"

mkdir -p "${FINAL_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${SAMTOOLS_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -c 'import sys; print("Python:", sys.version.replace("\n", " ")); raise SystemExit(0 if sys.version_info >= (3, 6) else "ERROR: Python >= 3.6 is required.")'

echo "=== Packaging frozen final genome ==="
echo "Source dir : ${SOURCE_DIR}"
echo "Final dir  : ${FINAL_DIR}"
echo "Prefix     : ${PREFIX}"
echo

samtools faidx "${FINAL_FASTA}"
samtools faidx "${CHR_FASTA}"
samtools faidx "${UNPLACED_FASTA}"
if [ -s "${COPY1_FASTA}" ]; then
    samtools faidx "${COPY1_FASTA}"
elif [ -s "${COPY1_FASTA_ALT}" ]; then
    samtools faidx "${COPY1_FASTA_ALT}"
fi
if [ -s "${COPY2_FASTA}" ]; then
    samtools faidx "${COPY2_FASTA}"
elif [ -s "${COPY2_FASTA_ALT}" ]; then
    samtools faidx "${COPY2_FASTA_ALT}"
fi

if [ "${INCLUDE_LARGE_EVIDENCE}" -eq 1 ]; then
    "${PYTHON_BIN}" "${PACKAGE_SCRIPT}" \
        --source-dir "${SOURCE_DIR}" \
        --final-dir "${FINAL_DIR}" \
        --prefix "${PREFIX}" \
        --include-large-evidence
else
    "${PYTHON_BIN}" "${PACKAGE_SCRIPT}" \
        --source-dir "${SOURCE_DIR}" \
        --final-dir "${FINAL_DIR}" \
        --prefix "${PREFIX}"
fi

echo
echo "=== Frozen final genome package ==="
echo "Final directory       : ${FINAL_DIR}"
echo "Primary genome FASTA  : ${FINAL_DIR}/genome/${PREFIX}.final.fa"
echo "Chromosome confidence : ${FINAL_DIR}/tables/${PREFIX}.chromosome_confidence.tsv"
echo "Final summary         : ${FINAL_DIR}/tables/${PREFIX}.final_genome_summary.tsv"
echo "Manifest              : ${FINAL_DIR}/MANIFEST.tsv"
echo "Checksums             : ${FINAL_DIR}/CHECKSUMS.md5"
echo "README                : ${FINAL_DIR}/README_FINAL_GENOME.md"
