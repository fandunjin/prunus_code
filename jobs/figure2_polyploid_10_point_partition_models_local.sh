#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Compare all shared/independent POInT partition models.
#
# Models:
#   1) PZI+PSP+PPA shared
#   2) PZI+PSP shared + PPA independent
#   3) PZI+PPA shared + PSP independent
#   4) PSP+PPA shared + PZI independent
#   5) PZI+PSP+PPA independent

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SHARED_ROOT="${SHARED_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_official_reference_exact_r13}"
INDEPENDENT_ROOT="${INDEPENDENT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_point_shared_vs_independent_r2/02_independent_runs}"
PILLARS="${PILLARS:-${SHARED_ROOT}/02_ancestral_order/point_pillars_optimized.txt}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/07_point_partition_models}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
POINT_BIN="${POINT_BIN:-${ROOT_DIR}/software/POInT/POInT}"
MODEL_FILE="${MODEL_FILE:-${ROOT_DIR}/software/POInT/models/WGD_bias_fix_model.txt}"
PREPARE="${ROOT_DIR}/scripts/fig2_prepare_point_pair_inputs.py"
SUMMARIZE="${ROOT_DIR}/scripts/fig2_summarize_point_partition_models.py"
PLOT="${ROOT_DIR}/scripts/plot_fig2_panel_d_point_partition_models.py"
POINT_TOL="${POINT_TOL:-1e-7}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${SHARED_ROOT}/04_summary/point_model_fit_all.tsv"
require_file "${PILLARS}"
require_file "${MODEL_FILE}"
require_file "${PREPARE}"
require_file "${SUMMARIZE}"
require_file "${PLOT}"
[[ -x "${POINT_BIN}" ]] || die "POInT binary unavailable: ${POINT_BIN}"
[[ -d "${INDEPENDENT_ROOT}/PZI" ]] || die "Independent POInT runs unavailable: ${INDEPENDENT_ROOT}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT"

mkdir -p "${OUT_ROOT}"
INPUT_ROOT="${OUT_ROOT}/01_pair_inputs"
RUN_ROOT="${OUT_ROOT}/02_pair_runs"
SUMMARY_ROOT="${OUT_ROOT}/03_summary"
PLOT_ROOT="${OUT_ROOT}/04_panel_d"
mkdir -p "${RUN_ROOT}" "${SUMMARY_ROOT}" "${PLOT_ROOT}"

"${PYTHON_BIN}" "${PREPARE}" \
    --inputs "${SHARED_ROOT}/01_inputs" \
    --optimized "${PILLARS}" \
    --outdir "${INPUT_ROOT}"

for pair in PZI_PSP PZI_PPA PSP_PPA; do
    pdir="${INPUT_ROOT}/${pair}"
    rdir="${RUN_ROOT}/${pair}"
    mkdir -p "${rdir}"
    left="${pair%_*}"
    right="${pair#*_}"
    if [[ -e "${rdir}/.external_running" && ! -s "${rdir}/species_tree.nex.out" ]]; then
        echo "[INFO] Waiting for externally started pair fit: ${pair}"
        while [[ ! -s "${rdir}/species_tree.nex.out" ]]; do
            sleep 10
        done
        rm -f "${rdir}/.external_running"
    fi
    if [[ -s "${rdir}/species_tree.nex.out" ]]; then
        echo "[OK] Reusing completed pair fit: ${pair}"
        continue
    fi
    (
        cd "${rdir}"
        "${POINT_BIN}" \
            -d:2 \
            -g:"${pdir}/${left}_POInT_geneorder.txt" \
            -g:"${pdir}/${right}_POInT_geneorder.txt" \
            -o:"${pdir}/pillars.txt" \
            -m:"${MODEL_FILE}" \
            -t:"${pdir}/species_tree.nex" \
            -p:"${rdir}/posterior_probs.txt" \
            -z -l:"${POINT_TOL}" \
            > "${rdir}/point.stdout" \
            2> "${rdir}/point.stderr"
    )
    require_file "${pdir}/species_tree.nex.out"
    cp "${pdir}/species_tree.nex.out" "${rdir}/species_tree.nex.out"
done

N_PILLARS=$(($(wc -l < "${PILLARS}") - 1))
"${PYTHON_BIN}" "${SUMMARIZE}" \
    --shared-fit "${SHARED_ROOT}/04_summary/point_model_fit_all.tsv" \
    --pair-root "${RUN_ROOT}" \
    --independent-root "${INDEPENDENT_ROOT}" \
    --pillars "${N_PILLARS}" \
    --out "${SUMMARY_ROOT}/five_model_comparison.tsv"

"${PYTHON_BIN}" "${PLOT}" \
    --summary "${SUMMARY_ROOT}/five_model_comparison.tsv" \
    --outdir "${PLOT_ROOT}"

cat > "${OUT_ROOT}/software_and_parameters.txt" <<EOF
method=official POInT v1.63
model=WGD_bias_fix_model.txt
pillars=${N_PILLARS}
shared_all=existing rooted three-taxon shared fit
pair_shared=two-taxon rooted shared fit plus independent singleton likelihood
all_independent=existing independent effective likelihood sum
pair_parameter_count=3 process + 3 rooted two-taxon branches + 4 independent singleton = 10
all_shared_parameter_count=3 process + 5 rooted three-taxon branches = 8
all_independent_parameter_count=3 species x (3 process + 1 effective terminal branch) = 12
bic_sample_unit=common ${N_PILLARS} pillars for all five models
EOF
printf 'Five POInT partition models completed.\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] Five-model POInT partition analysis completed: ${OUT_ROOT}"
