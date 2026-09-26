#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Fit POInT with physical A/B duplicate columns held fixed and render a
# separate tree-mapped panel.  The native POInT permutation model and outputs
# are left untouched by this job.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STRICT_ROOT="${STRICT_ROOT:-${ROOT_DIR}/0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT_DIR="${INPUT_DIR:-${STRICT_ROOT}/15_point_repartition_official_reference_exact/01_inputs}"
ORDER_DIR="${ORDER_DIR:-${STRICT_ROOT}/15_point_repartition_official_reference_exact/02_ancestral_order}"
OUT_ROOT="${OUT_ROOT:-${STRICT_ROOT}/25_point_fixed_ab_physical_tracks}"
RUN_DIR="${OUT_ROOT}/01_model_run"
NORMALIZED_DIR="${RUN_DIR}/normalized"
FIGURE_DIR="${OUT_ROOT}/02_tree_mapped"
POINT_SRC="${POINT_SRC:-${ROOT_DIR}/software/POInT/src}"
POINT_BIN="${POINT_BIN:-${ROOT_DIR}/software/POInT/POInT_fixedAB}"
FIXED_MODEL_TEMPLATE="${FIXED_MODEL_TEMPLATE:-${ROOT_DIR}/software/POInT/models/WGD_bias_fix_fixedAB_model.txt}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
POINT_TOL="${POINT_TOL:-1e-7}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required executable unavailable: $1"; }

require_command "${PYTHON_BIN}"
require_file "${INPUT_DIR}/PZI_POInT_geneorder.txt"
require_file "${INPUT_DIR}/PSP_POInT_geneorder.txt"
require_file "${INPUT_DIR}/PPA_POInT_geneorder.txt"
require_file "${ORDER_DIR}/point_pillars_optimized.txt"
require_file "${ORDER_DIR}/optimized_pillar_map.tsv"
require_file "${INPUT_DIR}/WGD_bias_fix_model.txt"
require_file "${FIXED_MODEL_TEMPLATE}"
require_file "${INPUT_DIR}/point_species_tree.nex"
require_file "${ROOT_DIR}/scripts/fig2_parse_point_results.py"
require_file "${ROOT_DIR}/scripts/fig2_make_point_optimized_model.py"
require_file "${ROOT_DIR}/scripts/plot_fig2_panel_d_point_tree_mapped.py"

PILLAR_COUNT=$(awk 'NR > 1 && NF {n++} END {print n + 0}' "${ORDER_DIR}/point_pillars_optimized.txt")
# 2026-09-14: the pillar count is data-dependent; the previous 13,204 value was
# tied to the pre-revision inputs.  Set EXPECTED_PILLARS to enforce a value.
EXPECTED_PILLARS="${EXPECTED_PILLARS:-}"
if [[ -n "${EXPECTED_PILLARS}" ]]; then
    [[ "${PILLAR_COUNT}" == "${EXPECTED_PILLARS}" ]] || die "Expected ${EXPECTED_PILLARS} optimized pillars, found ${PILLAR_COUNT}"
fi
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"

if [[ ! -x "${POINT_BIN}" ]]; then
    make -C "${POINT_SRC}" POInT_fixedAB
fi
[[ -x "${POINT_BIN}" ]] || die "Fixed A/B POInT binary unavailable: ${POINT_BIN}"

mkdir -p "${RUN_DIR}"
cp "${INPUT_DIR}/point_species_tree.nex" "${RUN_DIR}/point_species_tree.nex"

(
    cd "${RUN_DIR}"
    "${POINT_BIN}" \
        -d:2 \
        -g:"${INPUT_DIR}/PZI_POInT_geneorder.txt" \
        -g:"${INPUT_DIR}/PSP_POInT_geneorder.txt" \
        -g:"${INPUT_DIR}/PPA_POInT_geneorder.txt" \
        -o:"${ORDER_DIR}/point_pillars_optimized.txt" \
        -m:"${FIXED_MODEL_TEMPLATE}" \
        -t:"${RUN_DIR}/point_species_tree.nex" \
        -p:"${RUN_DIR}/posterior_probs.txt" \
        -z -l:"${POINT_TOL}" -fixedAB \
        > "${RUN_DIR}/point.stdout" \
        2> "${RUN_DIR}/point.stderr"
)
require_file "${RUN_DIR}/posterior_probs.txt"
require_file "${RUN_DIR}/point_species_tree.nex.out"

cp "${RUN_DIR}/point_species_tree.nex.out" "${RUN_DIR}/conditional_tree.nex"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_make_point_optimized_model.py" \
    --template "${FIXED_MODEL_TEMPLATE}" \
    --tree-output "${RUN_DIR}/point_species_tree.nex.out" \
    --out "${RUN_DIR}/optimized_model.txt" \
    --fixed-ab

(
    cd "${RUN_DIR}"
    "${POINT_BIN}" \
        -d:2 \
        -g:"${INPUT_DIR}/PZI_POInT_geneorder.txt" \
        -g:"${INPUT_DIR}/PSP_POInT_geneorder.txt" \
        -g:"${INPUT_DIR}/PPA_POInT_geneorder.txt" \
        -o:"${ORDER_DIR}/point_pillars_optimized.txt" \
        -m:"${RUN_DIR}/optimized_model.txt" \
        -t:"${RUN_DIR}/conditional_tree.nex" \
        -p:"${RUN_DIR}/posterior_probs.txt" \
        -c:"${RUN_DIR}/conditional_probs.txt" \
        -no_opt -z -l:"${POINT_TOL}" -fixedAB \
        > "${RUN_DIR}/conditional.stdout" \
        2> "${RUN_DIR}/conditional.stderr"
)
require_file "${RUN_DIR}/conditional_probs.txt"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_parse_point_results.py" \
    --posterior "${RUN_DIR}/posterior_probs.txt" \
    --conditional "${RUN_DIR}/conditional_probs.txt" \
    --stdout "${RUN_DIR}/conditional.stdout" \
    --tree-output "${RUN_DIR}/point_species_tree.nex.out" \
    --model-name "reference_WGD_bias_fix_fixedAB" \
    --model-parameters 2 \
    --pillar-map "${ORDER_DIR}/optimized_pillar_map.tsv" \
    --outdir "${NORMALIZED_DIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_panel_d_point_tree_mapped.py" \
    --fit "${NORMALIZED_DIR}/point_model_fit.tsv" \
    --transition "${NORMALIZED_DIR}/point_branch_transition_summary.tsv" \
    --tree "${RUN_DIR}/point_species_tree.nex.out" \
    --model "${RUN_DIR}/optimized_model.txt" \
    --fixed-ab \
    --outdir "${FIGURE_DIR}"

{
    echo "method=POInT v1.63 with optional fixed A/B identity permutation"
    echo "mode=physical A/B duplicate columns; permutation state count fixed at 1"
    echo "physical_tracking=direct duplicate-column tracks with within-column adjacency only"
    echo "model_parameters=FixRate,Copy1Bias; SwitchProb omitted because permutation is fixed"
    echo "topology=PPA,(PSP,PZI)"
    echo "pillars=${PILLAR_COUNT}"
    echo "model=reference_WGD_bias_fix_fixedAB"
    echo "model_template=${FIXED_MODEL_TEMPLATE}"
    echo "run_dir=${RUN_DIR}"
    echo "figure_dir=${FIGURE_DIR}"
    echo "standard_native_outputs=not overwritten"
} > "${OUT_ROOT}/software_and_parameters.txt"
printf 'fixed A/B POInT calculation and tree-mapped visualization completed\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] Fixed A/B POInT analysis completed: ${OUT_ROOT}"
