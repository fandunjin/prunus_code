#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/copy2 assemblies are kept as their existing separate inputs below.
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
GENOME_DIR="${HAP_RENAME_DIR}/genome"
COPY_DIR="${COPY_DIR:-${ROOT_DIR}/final_assembly/haplotype_provenance}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/haploid/final_haplotypes/gc_merqury}"
STATS_SCRIPT="${ROOT_DIR}/scripts/haplotype_gc_heterozygosity.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MERQURY_ENV="${MERQURY_ENV:-merqury}"
MERQURY_SH="${MERQURY_SH:-}"
THREADS="${THREADS:-32}"
KMER="${KMER:-21}"

# Preferred input: an existing meryl database built from accurate WGS reads.
# If MERYL_DB is empty, set READS to space-separated PCR-free Illumina or HiFi files.
# Standard ONT and Hi-C reads must not be used for publication-grade Merqury metrics.
MERYL_DB="${MERYL_DB:-}"
READS="${READS:-}"

A1_FASTA="${A1_FASTA:-${GENOME_DIR}/pzi_hapA1.fa}"
A2_FASTA="${A2_FASTA:-${GENOME_DIR}/pzi_hapA2.fa}"
B1_FASTA="${B1_FASTA:-${GENOME_DIR}/pzi_hapB1.fa}"
B2_FASTA="${B2_FASTA:-${GENOME_DIR}/pzi_hapB2.fa}"
COPY1_FASTA="${COPY1_FASTA:-${COPY_DIR}/pzi_2n32.copy1.hap_provenance.chr16.fa}"
COPY2_FASTA="${COPY2_FASTA:-${COPY_DIR}/pzi_2n32.copy2.hap_provenance.chr16.fa}"

require_file() {
    [[ -s "$1" ]] || { echo "ERROR: missing or empty $2: $1" >&2; exit 1; }
}

require_path() {
    [[ -e "$1" ]] || { echo "ERROR: missing $2: $1" >&2; exit 1; }
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: required command is unavailable in ${MERQURY_ENV}: $1" >&2
        exit 1
    }
}

require_file "$CONDA_SH" "Conda initialization script"
require_file "$STATS_SCRIPT" "statistics helper"
require_file "$A1_FASTA" "A1 FASTA"
require_file "$A2_FASTA" "A2 FASTA"
require_file "$B1_FASTA" "B1 FASTA"
require_file "$B2_FASTA" "B2 FASTA"
require_file "$COPY1_FASTA" "copy1 FASTA"
require_file "$COPY2_FASTA" "copy2 FASTA"

source "$CONDA_SH"
conda activate "$MERQURY_ENV"
require_command python
require_command meryl
require_command Rscript

if [[ -z "$MERQURY_SH" ]]; then
    MERQURY_SH="$(command -v merqury.sh || true)"
fi
if [[ -z "$MERQURY_SH" && -n "${MERQURY:-}" && -s "${MERQURY}/merqury.sh" ]]; then
    MERQURY_SH="${MERQURY}/merqury.sh"
fi
require_file "$MERQURY_SH" "Merqury driver script"

export OMP_NUM_THREADS="$THREADS"
mkdir -p "$OUT_DIR" "$OUT_DIR/meryl" "$OUT_DIR/single" "$OUT_DIR/pairs" "$OUT_DIR/tables"

if [[ -z "$MERYL_DB" ]]; then
    for candidate in \
        "${ROOT_DIR}/QC/read.meryl" \
        "${ROOT_DIR}/final_assembly/qc/merqury/read.meryl" \
        "${ROOT_DIR}/merqury/read.meryl"; do
        if [[ -d "$candidate" ]]; then
            MERYL_DB="$candidate"
            break
        fi
    done
fi

if [[ -z "$MERYL_DB" ]]; then
    if [[ -z "$READS" ]]; then
        echo "ERROR: Merqury requires accurate PCR-free Illumina/HiFi reads." >&2
        echo "Set MERYL_DB=/path/read.k${KMER}.meryl or READS='/path/R1.fq.gz /path/R2.fq.gz'." >&2
        echo "Do not use the ONT assembly reads or Hi-C reads for publication-grade QV/completeness." >&2
        exit 1
    fi
    read -r -a READ_ARRAY <<< "$READS"
    for read_file in "${READ_ARRAY[@]}"; do
        require_file "$read_file" "accurate WGS read file"
    done
    MERYL_DB="$OUT_DIR/meryl/pzi_reads.k${KMER}.meryl"
    if [[ ! -d "$MERYL_DB" ]]; then
        echo "[INFO] Building meryl read database: $MERYL_DB"
        meryl k="$KMER" count output "$MERYL_DB" "${READ_ARRAY[@]}"
    fi
fi
require_path "$MERYL_DB" "meryl read database"

GC_TABLE="$OUT_DIR/tables/gc_summary.tsv"
SUMMARY="$OUT_DIR/tables/gc_merqury_summary.tsv"
python "$STATS_SCRIPT" gc \
    --fasta "A1=$A1_FASTA" --fasta "A2=$A2_FASTA" \
    --fasta "B1=$B1_FASTA" --fasta "B2=$B2_FASTA" \
    --output "$GC_TABLE"

LABELS=(A1 A2 B1 B2 copy1 copy2)
declare -A FASTAS=(
    [A1]="$A1_FASTA" [A2]="$A2_FASTA" [B1]="$B1_FASTA" [B2]="$B2_FASTA"
    [copy1]="$COPY1_FASTA" [copy2]="$COPY2_FASTA"
)
RUN_ARGS=()

echo "=== PZI GC + Merqury ==="
echo "Targets    : ${LABELS[*]}"
echo "Meryl DB   : $MERYL_DB"
echo "k-mer      : $KMER"
echo "Merqury    : $MERQURY_SH"
echo "Output     : $OUT_DIR"

for label in "${LABELS[@]}"; do
    run_dir="$OUT_DIR/single/$label"
    prefix="$run_dir/$label"
    mkdir -p "$run_dir"
    if [[ ! -s "$prefix.qv" || ! -s "$prefix.completeness.stats" ]]; then
        echo "[INFO] Merqury single assembly: $label"
        bash "$MERQURY_SH" "$MERYL_DB" "${FASTAS[$label]}" "$prefix"
    else
        echo "[INFO] Existing Merqury result; skip: $label"
    fi
    require_file "$prefix.qv" "$label Merqury QV"
    require_file "$prefix.completeness.stats" "$label Merqury completeness"
    RUN_ARGS+=(--run "$label=$prefix")
done

run_pair() {
    local label="$1"
    local first="$2"
    local second="$3"
    local pair_dir="$OUT_DIR/pairs/$label"
    local prefix="$pair_dir/$label"
    mkdir -p "$pair_dir"
    if [[ ! -s "$prefix.spectra-asm.hist" ]]; then
        echo "[INFO] Merqury assembly pair: $label"
        bash "$MERQURY_SH" "$MERYL_DB" "$first" "$second" "$prefix"
    else
        echo "[INFO] Existing Merqury pair result; skip: $label"
    fi
}

run_pair A1_vs_A2 "$A1_FASTA" "$A2_FASTA"
run_pair B1_vs_B2 "$B1_FASTA" "$B2_FASTA"
run_pair copy1_vs_copy2 "$COPY1_FASTA" "$COPY2_FASTA"

python "$STATS_SCRIPT" merqury-summary \
    --gc "$GC_TABLE" "${RUN_ARGS[@]}" --output "$SUMMARY"

echo
echo "=== Results ==="
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo
echo "GC table       : $GC_TABLE"
echo "Merqury table  : $SUMMARY"
echo "Single runs    : $OUT_DIR/single"
echo "Pair spectra   : $OUT_DIR/pairs"
echo
echo "NOTE: Merqury error rate is consensus k-mer error, not biological heterozygosity."
echo "Use spectra-cn and spectra-asm from pair runs to assess residual duplication and haplotype separation."
