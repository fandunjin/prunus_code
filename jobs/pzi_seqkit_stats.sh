#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
#
# Optional overrides:
#   ROOT_DIR=/path/to/prunus RAW_DIR=/path/to/rawdata \
#   FINAL_FASTA=/path/to/pzi_2n32.final.fa OUT_DIR=/path/to/output \
#   jsub < jobs/pzi_seqkit_stats.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
SEQKIT_ENV="${SEQKIT_ENV:-seqkit}"
THREADS="${THREADS:-8}"

# Resolve both layouts used in this project. The first branch matches the
# project-root layout; the second matches ../rawdata/prunus/rawdata.
RAW_DIR="${RAW_DIR:-}"
if [[ -z "${RAW_DIR}" ]]; then
    if [[ -d "${ROOT_DIR}/rawdata/HD-14-2DNA" ]]; then
        RAW_DIR="${ROOT_DIR}/rawdata"
    elif [[ -d "${ROOT_DIR}/../rawdata/prunus/rawdata/HD-14-2DNA" ]]; then
        RAW_DIR="${ROOT_DIR}/../rawdata/prunus/rawdata"
    else
        RAW_DIR="${ROOT_DIR}/rawdata"
    fi
fi

ONT="${ONT:-${RAW_DIR}/HD-14-2DNA/HD-14-2DNA.pass.fastq.gz}"
HIC_R1="${HIC_R1:-${RAW_DIR}/HD-14-2DNA/H260100013A_1.fq.gz}"
HIC_R2="${HIC_R2:-${RAW_DIR}/HD-14-2DNA/H260100013A_2.fq.gz}"
ISOSEQ="${ISOSEQ:-${RAW_DIR}/HD-14-2PB/m84072_260707_083113_s2.fastq.gz}"
RNA_DIR="${RNA_DIR:-${RAW_DIR}/trans}"
FINAL_FASTA="${FINAL_FASTA:-${ROOT_DIR}/final_assembly/pzi_2n32.final.fa}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/analysis/pzi_seqkit}"

GENOME_STATS="${OUT_DIR}/pzi_genome_and_isoseq_seqkit.tsv"
RNA_STATS="${OUT_DIR}/pzi_rnaseq_seqkit.tsv"
ASSEMBLY_STATS="${OUT_DIR}/pzi_final_assembly_seqkit.tsv"
SUMMARY="${OUT_DIR}/pzi_seqkit_summary.tsv"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${ONT}" "ONT reads"
require_file "${HIC_R1}" "Hi-C R1 reads"
require_file "${HIC_R2}" "Hi-C R2 reads"
require_file "${ISOSEQ}" "PacBio Iso-Seq reads"
require_file "${FINAL_FASTA}" "final PZI assembly"
[[ -d "${RNA_DIR}" ]] || die "RNA-seq directory does not exist: ${RNA_DIR}"

shopt -s nullglob
RNA=("${RNA_DIR}"/*/*_1.fq.gz "${RNA_DIR}"/*/*_2.fq.gz)
(( ${#RNA[@]} > 0 )) || die "No Illumina RNA-seq FASTQ files found below: ${RNA_DIR}"
shopt -u nullglob

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${SEQKIT_ENV}"
command -v seqkit >/dev/null 2>&1 || die "seqkit is unavailable in environment: ${SEQKIT_ENV}"

echo "=== PZI seqkit statistics ==="
echo "Project root : ${ROOT_DIR}"
echo "Raw data     : ${RAW_DIR}"
echo "RNA files    : ${#RNA[@]}"
echo "Output       : ${OUT_DIR}"
echo

seqkit stats -a -T -j "${THREADS}" \
    "${ONT}" "${HIC_R1}" "${HIC_R2}" "${ISOSEQ}" \
    | tee "${GENOME_STATS}"

seqkit stats -a -T -j "${THREADS}" "${RNA[@]}" \
    | tee "${RNA_STATS}"

seqkit stats -a -T -j "${THREADS}" "${FINAL_FASTA}" \
    | tee "${ASSEMBLY_STATS}"

# Summarize the fields needed for reporting. Input paths are passed explicitly
# so that paths containing spaces or different absolute prefixes remain safe.
awk -F '\t' \
    -v ont="${ONT}" \
    -v hic_r1="${HIC_R1}" \
    -v hic_r2="${HIC_R2}" \
    -v isoseq="${ISOSEQ}" \
    -v assembly="${FINAL_FASTA}" \
    '
    function clean(value) {
        gsub(",", "", value)
        return value + 0
    }
    FNR == 1 { next }
    {
        n = clean($4)
        bp = clean($5)
        file = $1
        if (file == ont) {
            ont_reads = n; ont_bp = bp; ont_n50 = clean($13)
        } else if (file == hic_r1) {
            hic_r1_reads = n; hic_r1_bp = bp
        } else if (file == hic_r2) {
            hic_r2_reads = n; hic_r2_bp = bp
        } else if (file == isoseq) {
            isoseq_reads = n; isoseq_bp = bp; isoseq_n50 = clean($13)
        } else if (file == assembly) {
            assembly_seqs = n; assembly_bp = bp; assembly_n50 = clean($13)
        }
    }
    END {
        printf "metric\tvalue\n"
        printf "ONT_reads\t%.0f\n", ont_reads
        printf "ONT_total_bases\t%.0f\n", ont_bp
        printf "ONT_total_Gb\t%.3f\n", ont_bp / 1e9
        printf "HiC_R1_reads\t%.0f\n", hic_r1_reads
        printf "HiC_R2_reads\t%.0f\n", hic_r2_reads
        printf "HiC_read_pairs\t%.0f\n", hic_r1_reads
        printf "HiC_total_bases\t%.0f\n", hic_r1_bp + hic_r2_bp
        printf "HiC_total_Gb\t%.3f\n", (hic_r1_bp + hic_r2_bp) / 1e9
        printf "IsoSeq_reads\t%.0f\n", isoseq_reads
        printf "IsoSeq_total_bases\t%.0f\n", isoseq_bp
        printf "IsoSeq_total_Gb\t%.3f\n", isoseq_bp / 1e9
        printf "IsoSeq_read_N50\t%.0f\n", isoseq_n50
        printf "Final_assembly_sequences\t%.0f\n", assembly_seqs
        printf "Final_assembly_bases\t%.0f\n", assembly_bp
        printf "Final_assembly_Gb\t%.3f\n", assembly_bp / 1e9
        printf "Final_assembly_N50\t%.0f\n", assembly_n50
    }' "${GENOME_STATS}" "${ASSEMBLY_STATS}" > "${SUMMARY}.genome.tmp"

awk -F '\t' '
    FNR == 1 { next }
    {
        gsub(",", "", $4)
        gsub(",", "", $5)
        total_bp += $5 + 0
        if ($1 ~ /_1\.fq\.gz$/) {
            read_pairs += $4 + 0
        }
    }
    END {
        printf "RNA_total_bases\t%.0f\n", total_bp
        printf "RNA_total_Gb\t%.3f\n", total_bp / 1e9
        printf "RNA_read_pairs\t%.0f\n", read_pairs
    }' "${RNA_STATS}" >> "${SUMMARY}.genome.tmp"

{
    printf "metric\tvalue\n"
    tail -n +2 "${SUMMARY}.genome.tmp"
} > "${SUMMARY}"
rm -f "${SUMMARY}.genome.tmp"

echo
echo "Summary written to: ${SUMMARY}"
