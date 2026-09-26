#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/downstream_analysis/karyotype"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RIBBON_RUN_DIR:-${OUT_ROOT}/tree_pairwise_ribbon_chr3_11_12_13_15_${RUN_ID}}"

echo "=== JCVI tree-ordered pairwise ribbon ==="
echo "Highlighted pzi chromosomes: 3, 11, 12, 13, 15"
echo "Output: ${RUN_DIR}"

# Reuse the existing workflow, anchors, chromosome ordering and filtering.
exec env \
    HIGHLIGHT_ORDINALS="3 11 12 13 15" \
    HIGHLIGHT_COLORS="#0072B2 #D55E00 #009E73 #CC79A7 #E69F00" \
    RIBBON_RUN_DIR="${RUN_DIR}" \
    bash "${ROOT_DIR}/jobs/jcvi_tree_pairwise_ribbon.sh"
