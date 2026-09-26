#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/downstream_analysis/karyotype"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RIBBON_RUN_DIR:-${OUT_ROOT}/tree_pairwise_ribbon_chr1_chr2_${RUN_ID}}"

echo "=== JCVI tree-ordered pairwise ribbon: highlight pzi chromosomes 1 and 2 ==="
echo "Output: ${RUN_DIR}"

# Reuse the existing JCVI workflow and only change the highlighted pzi chromosomes.
exec env \
    HIGHLIGHT_ORDINALS="1 2" \
    HIGHLIGHT_COLORS="#56B4E9 #F0E442" \
    RIBBON_RUN_DIR="${RUN_DIR}" \
    bash "${ROOT_DIR}/jobs/jcvi_tree_pairwise_ribbon.sh"
