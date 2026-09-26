#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Refit the formal three-taxon POInT model and write native conditional
# branch-transition probabilities after the local POInT v1.63 memory fixes.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
POINT_ROOT="${POINT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_point_conditional_fixed_r1}"
PARTITION_ROOT="${PARTITION_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/07_point_partition_models_r1}"
PANEL_ROOT="${PANEL_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/10_point_panel_d_integrated_conditional}"

[[ ! -e "${POINT_ROOT}" ]] || {
    echo "[ERROR] Output directory already exists: ${POINT_ROOT}; set POINT_ROOT for a fresh run" >&2
    exit 2
}
[[ ! -e "${PANEL_ROOT}" ]] || {
    echo "[ERROR] Output directory already exists: ${PANEL_ROOT}; set PANEL_ROOT for a fresh run" >&2
    exit 2
}

POINT_CONDITIONAL=1 OUT_ROOT="${POINT_ROOT}" \
    bash "${ROOT_DIR}/jobs/figure2_polyploid_08_point_official_local.sh"

for model in symmetric reference_WGD_bias_fix; do
    conditional="${POINT_ROOT}/03_model_runs/${model}/conditional_probs.txt"
    [[ -s "${conditional}" ]] || {
        echo "[ERROR] Native POInT -c did not write ${conditional}" >&2
        exit 2
    }
done

POINT_ROOT="${POINT_ROOT}" PARTITION_ROOT="${PARTITION_ROOT}" OUT_ROOT="${PANEL_ROOT}" \
    bash "${ROOT_DIR}/jobs/figure2_polyploid_11_point_panel_d_integrated_local.sh"

echo "[OK] POInT conditional Figure 2d completed: ${PANEL_ROOT}"
