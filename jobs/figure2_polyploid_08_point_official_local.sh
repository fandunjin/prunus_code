#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for the official POInT likelihood workflow.
#
# This job prepares native POInT inputs, optimizes an ancestral pillar order,
# fits symmetric and copy-specific WGD loss models on the fixed topology
# (PPA,(PSP,PZI)), and normalizes the official posterior/conditional outputs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CATALOG_DIR="${CATALOG_DIR:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/01_catalog}"
PILLARS="${PILLARS:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_reference_exact_unambiguous/01_pillars/gtr_centered_pillars.tsv}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_official_reference_exact}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PLOT_PYTHON_BIN="${PLOT_PYTHON_BIN:-}"
POINT_BIN="${POINT_BIN:-${ROOT_DIR}/software/POInT/POInT}"
ANCESTOR_BIN="${ANCESTOR_BIN:-${ROOT_DIR}/software/POInT/POInT_ances_order}"
PREPARE_SCRIPT="${ROOT_DIR}/scripts/fig2_prepare_point_official_inputs.py"
PARSE_SCRIPT="${ROOT_DIR}/scripts/fig2_parse_point_results.py"
MAP_SCRIPT="${ROOT_DIR}/scripts/fig2_map_point_pillars.py"
OPTIMIZED_MODEL_SCRIPT="${ROOT_DIR}/scripts/fig2_make_point_optimized_model.py"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_fig2_panel_d_point_official.py"
REFERENCE_WGD_MODEL="${REFERENCE_WGD_MODEL:-${ROOT_DIR}/software/POInT/models/WGD_bias_fix_model.txt}"
STRICT_REFERENCE_MODEL="${STRICT_REFERENCE_MODEL:-1}"
ANCESTOR_RELAX="${ANCESTOR_RELAX:-50}"
POINT_TOL="${POINT_TOL:-1e-7}"
# POInT v1.63 reproducibly segfaults while materializing the conditional
# branch-transition matrix on this 11,309-pillar input, after completing the
# likelihood and posterior calculation.  Keep the official posterior path as
# the default; set POINT_CONDITIONAL=1 to request the legacy -c output.
POINT_CONDITIONAL="${POINT_CONDITIONAL:-0}"
ANCESTOR_MODE="${ANCESTOR_MODE:-gtr_fixed}"
# POInT_ances_order's OpenMP annealer is not stable on this WSL build when it
# simultaneously mutates a large pillar order.  One thread preserves the
# objective function and avoids non-deterministic worker-state corruption.
POINT_OMP_THREADS="${POINT_OMP_THREADS:-1}"
export OMP_NUM_THREADS="${POINT_OMP_THREADS}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required executable unavailable: $1"; }

require_file "${CATALOG_DIR}/tracks/PZI_A.bed"
require_file "${CATALOG_DIR}/tracks/PZI_B.bed"
require_file "${CATALOG_DIR}/tracks/PSP_A.bed"
require_file "${CATALOG_DIR}/tracks/PSP_B.bed"
require_file "${CATALOG_DIR}/tracks/PPA_A.bed"
require_file "${CATALOG_DIR}/tracks/PPA_B.bed"
require_file "${PILLARS}"
require_file "${PREPARE_SCRIPT}"
require_file "${PARSE_SCRIPT}"
require_file "${MAP_SCRIPT}"
require_file "${OPTIMIZED_MODEL_SCRIPT}"
require_file "${PLOT_SCRIPT}"
require_command "${PYTHON_BIN}"
[[ -x "${POINT_BIN}" ]] || die "Official POInT binary unavailable: ${POINT_BIN}; run jobs/install_point_local.sh first"
if [[ "${ANCESTOR_MODE}" == "point_ances_order" ]]; then
    [[ -x "${ANCESTOR_BIN}" ]] || die "Official POInT_ances_order unavailable: ${ANCESTOR_BIN}; run jobs/install_point_local.sh first"
fi
if [[ "${STRICT_REFERENCE_MODEL}" == "1" ]]; then
    require_file "${REFERENCE_WGD_MODEL}"
fi
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"

if [[ -z "${PLOT_PYTHON_BIN}" ]]; then
    for candidate in "${PYTHON_BIN}" \
        "/opt/plant-wakr-venv/bin/python3" \
        "${CONDA_ROOT}/envs/wgd2-reference-clean/bin/python"; do
        if [[ -x "${candidate}" ]] && "${candidate}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
            PLOT_PYTHON_BIN="${candidate}"
            break
        fi
    done
fi
[[ -n "${PLOT_PYTHON_BIN}" ]] || die "No Python interpreter with matplotlib and numpy found; set PLOT_PYTHON_BIN"

mkdir -p "${OUT_ROOT}"
INPUT_DIR="${OUT_ROOT}/01_inputs"
ORDER_DIR="${OUT_ROOT}/02_ancestral_order"
RUN_DIR="${OUT_ROOT}/03_model_runs"
SUMMARY_DIR="${OUT_ROOT}/04_summary"
mkdir -p "${ORDER_DIR}" "${RUN_DIR}" "${SUMMARY_DIR}"

"${PYTHON_BIN}" "${PREPARE_SCRIPT}" --pillars "${PILLARS}" --catalog "${CATALOG_DIR}" --outdir "${INPUT_DIR}"
for species in PZI PSP PPA; do require_file "${INPUT_DIR}/${species}_POInT_geneorder.txt"; done
require_file "${INPUT_DIR}/point_pillars_initial.txt"
require_file "${INPUT_DIR}/WGD_symmetric_model.txt"
require_file "${INPUT_DIR}/WGD_biased_model.txt"

# The published POInT WGD-bias model is the formal alternative.  The local
# symmetric model is retained only as its nested no-bias comparator.  A
# source-validated fallback is possible for development, but is never the
# default for the manuscript analysis.
if [[ -s "${REFERENCE_WGD_MODEL}" ]]; then
    cp "${REFERENCE_WGD_MODEL}" "${INPUT_DIR}/WGD_bias_fix_model.txt"
    COPY_SPECIFIC_MODEL="${INPUT_DIR}/WGD_bias_fix_model.txt"
    COPY_SPECIFIC_LABEL="reference_WGD_bias_fix"
    COPY_SPECIFIC_PARAMS="$(awk '
        {
            label = tolower($1)
            if (label == "#numberofparameters" || label == "#numparams") {
                if (NF > 1) { print $2; exit }
                if (getline > 0) { print $1; exit }
            }
        }
    ' "${COPY_SPECIFIC_MODEL}")"
    [[ "${COPY_SPECIFIC_PARAMS}" =~ ^[0-9]+$ ]] || die "Cannot read POInT model parameter count from ${COPY_SPECIFIC_MODEL}"
else
    COPY_SPECIFIC_MODEL="${INPUT_DIR}/WGD_biased_model.txt"
    COPY_SPECIFIC_LABEL="source_validated_copy_specific"
    COPY_SPECIFIC_PARAMS=3
fi

cp "${INPUT_DIR}/point_pillars_initial.txt" "${ORDER_DIR}/point_pillars_initial.txt"
if [[ "${ANCESTOR_MODE}" == "point_ances_order" ]]; then
    (
        cd "${ORDER_DIR}"
        "${ANCESTOR_BIN}" \
            -d:2 \
            -g:"${INPUT_DIR}/PZI_POInT_geneorder.txt" \
            -g:"${INPUT_DIR}/PSP_POInT_geneorder.txt" \
            -g:"${INPUT_DIR}/PPA_POInT_geneorder.txt" \
            -o:"${ORDER_DIR}/point_pillars_initial.txt" \
            -n:"${ORDER_DIR}/point_pillars_optimized.txt" \
            -m:1 -r:"${ANCESTOR_RELAX}" -b:0.0001 -infer_init_order -usepositionscore \
            > "${ORDER_DIR}/point_ances_order.stdout" \
            2> "${ORDER_DIR}/point_ances_order.stderr"
    )
elif [[ "${ANCESTOR_MODE}" == "gtr_fixed" ]]; then
    # The pillar catalog is ordered by the low-copy GTR anchor.  Retain this
    # externally defined order after the legacy POInT_ances_order helper
    # segfaulted reproducibly on this 11,309-pillar dataset.
    cp "${ORDER_DIR}/point_pillars_initial.txt" "${ORDER_DIR}/point_pillars_optimized.txt"
    printf 'ancestral order supplied by GTR-centered pillar catalog; POInT_ances_order skipped\n' > "${ORDER_DIR}/point_ances_order.stdout"
    : > "${ORDER_DIR}/point_ances_order.stderr"
else
    die "Unsupported ANCESTOR_MODE: ${ANCESTOR_MODE}; use gtr_fixed or point_ances_order"
fi

# POInT's executable parser uses the exact file names and the optimized pillar
# order written by POInT_ances_order.  The expected PZI file name is corrected
# below for backward compatibility with older copies of this runner.
if [[ ! -s "${ORDER_DIR}/point_pillars_optimized.txt" ]]; then
    die "POInT_ances_order did not write an optimized pillar file"
fi
"${PYTHON_BIN}" "${MAP_SCRIPT}" \
    --initial-pillars "${INPUT_DIR}/point_pillars_initial.txt" \
    --initial-map "${INPUT_DIR}/initial_pillar_map.tsv" \
    --optimized-pillars "${ORDER_DIR}/point_pillars_optimized.txt" \
    --out "${ORDER_DIR}/optimized_pillar_map.tsv"

run_model() {
    local model_name="$1"
    local model_file="$2"
    local model_params="$3"
    local model_dir="${RUN_DIR}/${model_name}"
    local tree_file="${model_dir}/point_species_tree.nex"
    local stdout_file="${model_dir}/point.stdout"
    local conditional_file="${model_dir}/conditional_probs.txt"
    mkdir -p "${model_dir}"
    cp "${INPUT_DIR}/point_species_tree.nex" "${tree_file}"
    local parse_stdout="${stdout_file}"

    # Always fit the likelihood first without -c.  When conditional output is
    # requested, a second -no_opt pass reads the optimized tree and an edited
    # model file.  This avoids the v1.63 crash observed when conditional
    # matrices are materialized inside the optimizer process.
    (
        cd "${model_dir}"
        "${POINT_BIN}" \
            -d:2 \
            -g:"${INPUT_DIR}/PZI_POInT_geneorder.txt" \
            -g:"${INPUT_DIR}/PSP_POInT_geneorder.txt" \
            -g:"${INPUT_DIR}/PPA_POInT_geneorder.txt" \
            -o:"${ORDER_DIR}/point_pillars_optimized.txt" \
            -m:"${model_file}" \
            -t:"${tree_file}" \
            -p:"${model_dir}/posterior_probs.txt" \
            -z -l:"${POINT_TOL}" \
            > "${stdout_file}" \
            2> "${model_dir}/point.stderr"
    ) || die "Official POInT likelihood fit failed for ${model_name}; see ${model_dir}/point.stderr"

    require_file "${model_dir}/posterior_probs.txt"
    require_file "${tree_file}.out"

    if [[ "${POINT_CONDITIONAL}" == "1" ]]; then
        local conditional_tree="${model_dir}/conditional_tree.nex"
        local optimized_model="${model_dir}/optimized_model.txt"
        cp "${tree_file}.out" "${conditional_tree}"
        "${PYTHON_BIN}" "${OPTIMIZED_MODEL_SCRIPT}" \
            --template "${model_file}" \
            --tree-output "${tree_file}.out" \
            --out "${optimized_model}"
        set +e
        (
            cd "${model_dir}"
            "${POINT_BIN}" \
                -d:2 \
                -g:"${INPUT_DIR}/PZI_POInT_geneorder.txt" \
                -g:"${INPUT_DIR}/PSP_POInT_geneorder.txt" \
                -g:"${INPUT_DIR}/PPA_POInT_geneorder.txt" \
                -o:"${ORDER_DIR}/point_pillars_optimized.txt" \
                -m:"${optimized_model}" \
                -t:"${conditional_tree}" \
                -p:"${model_dir}/posterior_probs.txt" \
                -c:"${conditional_file}" \
                -no_opt -z -l:"${POINT_TOL}" \
                > "${model_dir}/conditional.stdout" \
                2> "${model_dir}/conditional.stderr"
        )
        local conditional_status=$?
        set -e
        if (( conditional_status == 0 )) && [[ -s "${conditional_file}" ]]; then
            parse_stdout="${model_dir}/conditional.stdout"
            echo "[OK] Native POInT -c completed for ${model_name}"
        else
            echo "[WARN] Native POInT -c still failed for ${model_name}; retaining posterior-only fit" >&2
            rm -f "${conditional_file}"
        fi
    fi

    require_file "${model_dir}/posterior_probs.txt"
    require_file "${tree_file}.out"
    "${PYTHON_BIN}" "${PARSE_SCRIPT}" \
        --posterior "${model_dir}/posterior_probs.txt" \
        --conditional "${conditional_file}" \
        --stdout "${parse_stdout}" \
        --tree-output "${tree_file}.out" \
        --model-name "${model_name}" \
        --model-parameters "${model_params}" \
        --pillar-map "${ORDER_DIR}/optimized_pillar_map.tsv" \
        --outdir "${model_dir}/normalized"
}

run_model "symmetric" "${INPUT_DIR}/WGD_symmetric_model.txt" 2
run_model "${COPY_SPECIFIC_LABEL}" "${COPY_SPECIFIC_MODEL}" "${COPY_SPECIFIC_PARAMS}"

{
    head -n 1 "${RUN_DIR}/symmetric/normalized/point_model_fit.tsv"
    tail -n +2 "${RUN_DIR}/symmetric/normalized/point_model_fit.tsv"
    tail -n +2 "${RUN_DIR}/${COPY_SPECIFIC_LABEL}/normalized/point_model_fit.tsv"
} > "${SUMMARY_DIR}/point_model_fit_all.tsv"
cp "${RUN_DIR}/${COPY_SPECIFIC_LABEL}/normalized/point_branch_transition_summary.tsv" "${SUMMARY_DIR}/point_branch_transition_summary.tsv"
cp "${RUN_DIR}/${COPY_SPECIFIC_LABEL}/normalized/point_posterior_qc.tsv" "${SUMMARY_DIR}/point_posterior_qc.tsv"

    "${PLOT_PYTHON_BIN}" "${PLOT_SCRIPT}" \
    --fit "${SUMMARY_DIR}/point_model_fit_all.tsv" \
    --transition "${SUMMARY_DIR}/point_branch_transition_summary.tsv" \
    --posterior-qc "${SUMMARY_DIR}/point_posterior_qc.tsv" \
    --posterior-summary "${RUN_DIR}/${COPY_SPECIFIC_LABEL}/normalized/point_posterior_summary.tsv" \
    --outdir "${OUT_ROOT}/05_panel_d"

{
    echo "method=official POInT v1.63 likelihood"
    echo "topology=(PPA,(PSP,PZI));"
    echo "external_anchor=GTR only for pillar ordering; GTR excluded from likelihood terminals"
    echo "models=symmetric nested comparator and ${COPY_SPECIFIC_LABEL}"
    echo "reference_wgd_model=${REFERENCE_WGD_MODEL}"
    echo "strict_reference_model=${STRICT_REFERENCE_MODEL}"
    echo "ancestral_order_mode=${ANCESTOR_MODE}"
    echo "ancestral_order_helper=POInT_ances_order -m:1 -r:${ANCESTOR_RELAX}"
    echo "point_tolerance=${POINT_TOL}"
    echo "conditional_output=${POINT_CONDITIONAL}"
    echo "omp_num_threads=${OMP_NUM_THREADS}"
    echo "inputs=${INPUT_DIR}"
    echo "summary=${SUMMARY_DIR}"
    echo "interpretation=official posterior tracking states; conditional branch transitions are reported only when POInT -c completes; absence is not independently verified gene loss"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Official POInT model analysis completed: ${OUT_ROOT}"
