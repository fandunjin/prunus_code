#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/0903/haploid_repartition/inputs/psp/genome}"
EDTA_ROOT="${EDTA_ROOT:-${ROOT_DIR}/data/psp/haploid/annotation/repeats}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

case "$LABEL" in
    hapA2|hapB1|hapB2) ;;
    *) echo "ERROR: LABEL must be hapA2, hapB1 or hapB2" >&2; exit 2 ;;
esac

GENOME="${INPUT_ROOT}/psp_${LABEL}.fa"
OUT_DIR="${EDTA_ROOT}/${LABEL}_edta"
WORK_FASTA="${OUT_DIR}/psp_${LABEL}.fa"

[[ -s "$GENOME" ]] || { echo "ERROR: missing genome: $GENOME" >&2; exit 2; }
[[ -s "$EDTA_CDS" ]] || { echo "ERROR: missing CDS: $EDTA_CDS" >&2; exit 2; }
[[ -s "$CONDA_SH" ]] || { echo "ERROR: missing conda init: $CONDA_SH" >&2; exit 2; }
mkdir -p "$OUT_DIR" "${ROOT_DIR}/logs"

if [[ ! -s "$WORK_FASTA" ]]; then
    cp -p "$GENOME" "$WORK_FASTA"
elif ! cmp -s "$GENOME" "$WORK_FASTA"; then
    echo "ERROR: work FASTA differs from input: $WORK_FASTA" >&2
    exit 2
fi

source "$CONDA_SH"
conda activate "$EDTA_ENV"
cd "$OUT_DIR"

echo "=== PSP independent EDTA completion ==="
echo "Label  : $LABEL"
echo "Genome : $WORK_FASTA"
echo "Output : $OUT_DIR"
echo "CDS    : $EDTA_CDS"
echo "CPUs   : $THREADS"

EDTA.pl \
    --genome "$WORK_FASTA" \
    --species others \
    --step all \
    --anno 1 \
    --sensitive 1 \
    --force "$EDTA_FORCE" \
    --overwrite "$EDTA_OVERWRITE" \
    --threads "$THREADS" \
    --cds "$EDTA_CDS"

for output in \
    "${WORK_FASTA}.mod.EDTA.TEanno.gff3" \
    "${WORK_FASTA}.mod.EDTA.TEanno.sum" \
    "${WORK_FASTA}.mod.EDTA.TElib.fa" \
    "${WORK_FASTA}.mod.EDTA.intact.gff3" \
    "${WORK_FASTA}.mod.EDTA.intact.fa"; do
    [[ -s "$output" ]] || { echo "ERROR: missing EDTA output: $output" >&2; exit 1; }
done

echo "[OK] Complete independent EDTA outputs: $OUT_DIR"
