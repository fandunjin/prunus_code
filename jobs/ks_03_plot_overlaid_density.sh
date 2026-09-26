#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
KAKS_ENV="${KAKS_ENV:-akr}"
RESULT_DIR="${RESULT_DIR:-${KS_ROOT}/kaks/collinear_all_pairs_YN}"
INPUT_TSV="${INPUT_TSV:-${RESULT_DIR}/summary/all_merged_4dtv_kaks.tsv}"
OUT_DIR="${OUT_DIR:-${RESULT_DIR}/figures/formal_density}"
MIN_KS="${MIN_KS:-0}"
MAX_KS="${MAX_KS:-3}"
KDE_BANDWIDTH="${KDE_BANDWIDTH:-0.02}"
ALLOW_LEGACY_INPUT="${ALLOW_LEGACY_INPUT:-0}"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_ks_overlaid_density.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

mkdir -p "${KS_ROOT}/logs" "${OUT_DIR}"
rm -f "${OUT_DIR}/COMPLETE"
require_file "${INPUT_TSV}"
require_file "${PLOT_SCRIPT}"

set +u
source "${CONDA_SH}"
conda activate "${KAKS_ENV}"
set -u

for module in numpy matplotlib; do
    python -c "import ${module}" >/dev/null 2>&1 || die "Python module unavailable in ${KAKS_ENV}: ${module}"
done

extra_args=()
if [[ "${ALLOW_LEGACY_INPUT}" == "1" ]]; then
    extra_args+=(--allow-legacy-input)
fi

python "${PLOT_SCRIPT}" \
    --input "${INPUT_TSV}" \
    --outdir "${OUT_DIR}" \
    --presets focal sensitivity cross_species diploid_references \
    --min-ks "${MIN_KS}" \
    --max-ks "${MAX_KS}" \
    --bandwidth "${KDE_BANDWIDTH}" \
    --bandwidth-sensitivity 0.02 0.04 0.06 \
    "${extra_args[@]}"

for preset in focal sensitivity cross_species diploid_references; do
    require_file "${OUT_DIR}/ks_density_${preset}.svg"
    require_file "${OUT_DIR}/ks_density_${preset}.pdf"
    require_file "${OUT_DIR}/ks_density_${preset}.png"
    require_file "${OUT_DIR}/ks_density_${preset}.tiff"
    require_file "${OUT_DIR}/ks_density_${preset}.source_data.tsv"
done

for extension in svg pdf png tiff source_data.tsv summary.tsv; do
    require_file "${OUT_DIR}/ks_density_bandwidth_sensitivity.${extension}"
    require_file "${OUT_DIR}/ks_density_focal_with_references.${extension}"
done
require_file "${OUT_DIR}/formal_figure_manifest.tsv"

{
    echo "python=$(python --version 2>&1)"
    echo "numpy=$(python -c 'import numpy; print(numpy.__version__)')"
    echo "matplotlib=$(python -c 'import matplotlib; print(matplotlib.__version__)')"
    echo "input=$(basename "${INPUT_TSV}")"
    echo "ks_filter=${MIN_KS}<Ks<${MAX_KS}"
    echo "kde_bandwidth=${KDE_BANDWIDTH}"
    echo "kde_bandwidth_sensitivity=0.02,0.04,0.06"
    echo "allow_legacy_input=${ALLOW_LEGACY_INPUT}"
} > "${OUT_DIR}/software_and_parameters.txt"

touch "${OUT_DIR}/COMPLETE"
echo "[OK] Ks density plotting completed: ${OUT_DIR}"
