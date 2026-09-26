#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
STATS_DIR="${STATS_DIR:-${PAN_ROOT}/followup_v2}"
FIG_DIR="${FIG_DIR:-${PAN_ROOT}/figures_followup_v2}"
# Use the original unmasked PZI assemblies for monomer extraction.  Do not
# silently read a masked FASTA because N bases would create artificial logo
# gaps.  GENOME_ROOT may be overridden with an equivalent FASTA root.
GENOME_ROOT="${GENOME_ROOT:-${ROOT_DIR}/final_haplotype_genomes}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MAFFT_ENV="${MAFFT_ENV:-orthofinder}"
MAFFT_BIN="${MAFFT_BIN:-}"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_prunus_pan_te_followup.py"
TRUE_LOGO_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_trf_true_logo.py"
ARRAYS="${ARRAYS:-${PAN_ROOT}/figure1_tandem/score50/tandem_top_arrays.tsv}"

[ -s "${PLOT_SCRIPT}" ] || { echo "ERROR: missing ${PLOT_SCRIPT}" >&2; exit 1; }
[ -s "${TRUE_LOGO_SCRIPT}" ] || { echo "ERROR: missing ${TRUE_LOGO_SCRIPT}" >&2; exit 1; }
[ -s "${STATS_DIR}/te_trf_1mb_density.tsv" ] || { echo "ERROR: missing ${STATS_DIR}/te_trf_1mb_density.tsv" >&2; exit 1; }
[ -s "${ARRAYS}" ] || { echo "ERROR: missing long-array table ${ARRAYS}" >&2; exit 1; }
source "${CONDA_SH}"
conda activate haphic
export MPLBACKEND=Agg
PYTHON_BIN="${PYTHON_BIN:-$(command -v python)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
[ -x "${PYTHON_BIN}" ] || { echo "ERROR: haphic Python was not found" >&2; exit 1; }
[ -x "${SAMTOOLS_BIN}" ] || { echo "ERROR: samtools is required to index genome FASTAs" >&2; exit 1; }

# MAFFT is maintained in a separate environment on this cluster.  Keep the
# plotting/runtime Python from haphic, but resolve MAFFT independently.
if [ -z "${MAFFT_BIN}" ] && command -v mafft >/dev/null 2>&1; then
    MAFFT_BIN="$(command -v mafft)"
fi
if [ -z "${MAFFT_BIN}" ]; then
    conda activate "${MAFFT_ENV}"
    MAFFT_BIN="$(command -v mafft || true)"
fi
[ -x "${MAFFT_BIN}" ] || {
    echo "ERROR: MAFFT was not found. Set MAFFT_BIN=/path/to/mafft or MAFFT_ENV=<conda-env>." >&2
    echo "Available environments can be listed with: conda env list" >&2
    exit 1
}
echo "MAFFT: ${MAFFT_BIN}"
echo "Python: ${PYTHON_BIN}"
echo "samtools: ${SAMTOOLS_BIN}"
THREADS="${THREADS:-8}"
declare -A GENOME_FASTA=(
    [pziA1]="${GENOME_ROOT}/pzi_hapA1.fa"
    [pziA2]="${GENOME_ROOT}/pzi_hapA2.fa"
    [pziB1]="${GENOME_ROOT}/pzi_hapB1.fa"
    [pziB2]="${GENOME_ROOT}/pzi_hapB2.fa"
)
for target in pziA1 pziA2 pziB1 pziB2; do
    fasta="${GENOME_FASTA[$target]}"
    [ -s "${fasta}" ] || { echo "ERROR: missing unmasked genome ${fasta}" >&2; exit 1; }
    [ -s "${fasta}.fai" ] || "${SAMTOOLS_BIN}" faidx "${fasta}"
done
"${PYTHON_BIN}" "${TRUE_LOGO_SCRIPT}" \
    --arrays "${ARRAYS}" \
    --motif-clusters "${STATS_DIR}/trf_motif_clusters.tsv" \
    --motif-members "${STATS_DIR}/trf_motif_members.tsv" \
    --genome-root "${GENOME_ROOT}" \
    --output-dir "${STATS_DIR}" \
    --mafft "${MAFFT_BIN}" \
    --clusters 4 \
    --max-monomers-per-array 40 \
    --max-monomers-per-cluster 400 \
    --align-all \
    --min-seed-identity 0.80 \
    --coordinate-system zero-based-half-open \
    --threads "${THREADS}"
"${PYTHON_BIN}" "${PLOT_SCRIPT}" --stats-dir "${STATS_DIR}" --output-dir "${FIG_DIR}"
echo "=== Follow-up figure replot completed ==="
