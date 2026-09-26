#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
KAKS_ENV="${KAKS_ENV:-akr}"
INPUT="${INPUT:-${KS_ROOT}/kaks/collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv}"
OUTDIR="${OUTDIR:-${KS_ROOT}/kaks/blockwise_median_YN}"
SCRIPT="${ROOT_DIR}/scripts/ks_blockwise_median_analysis.py"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
MIXTURE_BOOTSTRAPS="${MIXTURE_BOOTSTRAPS:-1000}"
MIN_PAIRS_PER_BLOCK="${MIN_PAIRS_PER_BLOCK:-5}"
SEED="${SEED:-20260815}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

mkdir -p "${KS_ROOT}/logs" "${OUTDIR}"
rm -f "${OUTDIR}/COMPLETE"
require_file "${INPUT}"
require_file "${SCRIPT}"

set +u
source "${CONDA_SH}"
conda activate "${KAKS_ENV}"
set -u

for module in numpy matplotlib; do
    python -c "import ${module}" >/dev/null 2>&1 || die "Python module unavailable in ${KAKS_ENV}: ${module}"
done

python "${SCRIPT}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    --min-ks 0 \
    --max-ks 3 \
    --min-pairs-per-block "${MIN_PAIRS_PER_BLOCK}" \
    --min-pairs-sensitivity 3 5 10 \
    --bootstrap-replicates "${BOOTSTRAPS}" \
    --mixture-bootstrap-replicates "${MIXTURE_BOOTSTRAPS}" \
    --max-components 3 \
    --mixture-initializations 30 \
    --bootstrap-initializations 8 \
    --min-component-blocks 5 \
    --min-component-fraction 0.10 \
    --kde-bandwidth 0.02 \
    --seed "${SEED}"

for output in \
    block_medians.all.tsv \
    block_medians.retained.tsv \
    blockwise_bootstrap_summary.tsv \
    blockwise_bootstrap_contrasts.tsv \
    minimum_block_size_sensitivity.tsv \
    mixture_model_selection.tsv \
    mixture_bootstrap_stability.tsv \
    blockwise_median_ks.pdf \
    blockwise_median_ks.svg \
    blockwise_median_ks.png \
    blockwise_median_ks.tiff \
    COMPLETE; do
    require_file "${OUTDIR}/${output}"
done

{
    echo "python=$(python --version 2>&1)"
    echo "numpy=$(python -c 'import numpy; print(numpy.__version__)')"
    echo "matplotlib=$(python -c 'import matplotlib; print(matplotlib.__version__)')"
    echo "input=${INPUT}"
    echo "bootstrap_unit=JCVI_block_id"
    echo "bootstrap_replicates=${BOOTSTRAPS}"
    echo "mixture_bootstrap_replicates=${MIXTURE_BOOTSTRAPS}"
    echo "min_pairs_per_block=${MIN_PAIRS_PER_BLOCK}"
    echo "mixture_scale=natural_log_Ks"
    echo "mixture_selection=BIC_primary_AIC_reported"
    echo "minimum_component_effective_blocks=max(5,10_percent_of_retained_blocks)"
    echo "seed=${SEED}"
} > "${OUTDIR}/software_versions.txt"

echo "[OK] Block-wise median Ks stage completed: ${OUTDIR}"
