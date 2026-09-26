#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Re-run the complete POInT fractionation workflow on the 0903 inputs.
# All outputs are isolated under 0903/downstream_analysis.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_ROOT="${DATA_ROOT:-${ROOT_DIR}/0903}"
INPUT_ROOT="${INPUT_ROOT:-${DATA_ROOT}/downstream_analysis/figure2_repartition_ks/inputs}"
BASE="${BASE:-${DATA_ROOT}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"

FRACTBIAS_OUT="${FRACTBIAS_OUT:-${BASE}/04_fractbias_8track_reference_exact}"
PILLAR_OUT="${PILLAR_OUT:-${BASE}/14_point_repartition_reference_exact/01_pillars}"
POINT_OUT="${POINT_OUT:-${BASE}/15_point_repartition_official_reference_exact}"
SHARED_OUT="${SHARED_OUT:-${BASE}/16_point_shared_independent_repartition}"
PARTITION_OUT="${PARTITION_OUT:-${BASE}/17_point_partition_models_repartition}"
INTEGRATED_OUT="${INTEGRATED_OUT:-${BASE}/19_point_panel_d_integrated_repartition}"
TREE_OUT="${TREE_OUT:-${BASE}/18_point_repartition_tree_mapped}"

if [[ -z "${PYTHON_BIN:-}" ]]; then
    if [[ -x "/opt/plant-wakr-venv/bin/python3" ]]; then
        PYTHON_BIN="/opt/plant-wakr-venv/bin/python3"
    else
        PYTHON_BIN="python3"
    fi
fi
if [[ -z "${PLOT_PYTHON_BIN:-}" && -x "/opt/plant-wakr-venv/bin/python3" ]] \
    && "/opt/plant-wakr-venv/bin/python3" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
    PLOT_PYTHON_BIN="/opt/plant-wakr-venv/bin/python3"
fi
THREADS="${THREADS:-12}"
POINT_CONDITIONAL="${POINT_CONDITIONAL:-1}"
ANCESTOR_MODE="${ANCESTOR_MODE:-gtr_fixed}"
export PYTHON_BIN PLOT_PYTHON_BIN

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_absent() { [[ ! -e "$1" ]] || die "Output already exists: $1; choose a new output root or remove it after checking"; }

require_file "${INPUT_ROOT}/gtr/gtr.bed"
require_file "${INPUT_ROOT}/pcam/pcam.bed"
require_file "${INPUT_ROOT}/ppe/ppe.bed"
require_file "${INPUT_ROOT}/pzi_hapA1/pzi_hapA1.bed"
require_file "${INPUT_ROOT}/psp_hapA1/psp_hapA1.bed"
require_file "${INPUT_ROOT}/ppa_hapA/ppa_hapA.bed"
require_file "${DATA_ROOT}/haploid_repartition/inputs/ANNOTATIONS_COMPLETE"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"

for output in "${FRACTBIAS_OUT}" "${PILLAR_OUT}" "${POINT_OUT}" "${SHARED_OUT}" "${PARTITION_OUT}" "${INTEGRATED_OUT}" "${TREE_OUT}"; do
    require_absent "${output}"
done

echo "[1/7] Formal GTR anchoring, Quota Align and FractBias"
INPUT_ROOT="${INPUT_ROOT}" \
OUT_ROOT="${FRACTBIAS_OUT}" \
THREADS="${THREADS}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_07_fractbias_reference_local.sh"

echo "[2/7] GTR-centered POInT pillar preparation"
mkdir -p "${PILLAR_OUT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_prepare_point_pillars.py" \
    --quota-anchors "${FRACTBIAS_OUT}/04_quota_align/normalized/quota_filtered_anchors.tsv" \
    --catalog "${FRACTBIAS_OUT}/01_catalog" \
    --outdir "${PILLAR_OUT}"

echo "[3/7] Official POInT shared-WGD likelihood and conditional states"
CATALOG_DIR="${FRACTBIAS_OUT}/01_catalog" \
PILLARS="${PILLAR_OUT}/gtr_centered_pillars.tsv" \
OUT_ROOT="${POINT_OUT}" \
POINT_CONDITIONAL="${POINT_CONDITIONAL}" \
ANCESTOR_MODE="${ANCESTOR_MODE}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_08_point_official_local.sh"

echo "[4/7] Shared-WGD versus independent-polyploidy likelihoods"
SHARED_ROOT="${POINT_OUT}" \
PILLARS="${POINT_OUT}/02_ancestral_order/point_pillars_optimized.txt" \
OUT_ROOT="${SHARED_OUT}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_09_point_shared_independent_local.sh"

echo "[5/7] Five shared/independent partition models and AIC/BIC"
SHARED_ROOT="${POINT_OUT}" \
INDEPENDENT_ROOT="${SHARED_OUT}/02_independent_runs" \
PILLARS="${POINT_OUT}/02_ancestral_order/point_pillars_optimized.txt" \
OUT_ROOT="${PARTITION_OUT}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_10_point_partition_models_local.sh"

echo "[6/7] Integrated POInT panel"
POINT_ROOT="${POINT_OUT}" \
PARTITION_ROOT="${PARTITION_OUT}" \
OUT_ROOT="${INTEGRATED_OUT}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_11_point_panel_d_integrated_local.sh"

echo "[7/7] Reference-style tree-mapped Figure 2d"
POINT_ROOT="${POINT_OUT}" \
PARTITION_ROOT="${PARTITION_OUT}" \
OUT_ROOT="${TREE_OUT}" \
bash "${ROOT_DIR}/jobs/figure2_polyploid_14_point_tree_mapped_existing_local.sh"

cat > "${BASE}/0903_point_run_manifest.txt" <<EOF
input_root=${INPUT_ROOT}
fractbias=${FRACTBIAS_OUT}
pillars=${PILLAR_OUT}
official_point=${POINT_OUT}
shared_independent=${SHARED_OUT}
partition_models=${PARTITION_OUT}
integrated_panel=${INTEGRATED_OUT}
tree_mapped_panel=${TREE_OUT}
point_conditional=${POINT_CONDITIONAL}
ancestor_mode=${ANCESTOR_MODE}
threads=${THREADS}
EOF

printf '0903 POInT workflow completed.\n' > "${BASE}/COMPLETE"
echo "[OK] Complete 0903 POInT workflow: ${BASE}"
echo "[OK] Final tree-mapped PDF: ${TREE_OUT}/figure2_panel_d_point_tree_mapped.pdf"
