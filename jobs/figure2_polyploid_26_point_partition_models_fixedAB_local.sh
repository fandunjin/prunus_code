#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# FixedAB variant of figure2_polyploid_10_point_partition_models_local.sh.
#
# Re-runs the five-model POInT partition comparison with the modified build that
# holds subgenome identity fixed at the physical A/B columns (POInT_fixedAB and
# WGD_bias_fix_fixedAB_model.txt).  Inputs are copied from the official-build
# comparison (17_point_partition_models_repartition) and the independent inputs
# from 16_point_shared_independent_repartition, so the two comparisons differ
# only in the build.
#
# Run from the repository root (WSL):
#   bash jobs/figure2_polyploid_26_point_partition_models_fixedAB_local.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STRICT_ROOT="${STRICT_ROOT:-${ROOT_DIR}/0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
PAIR_INPUTS="${PAIR_INPUTS:-${STRICT_ROOT}/17_point_partition_models_repartition/01_pair_inputs}"
IND_INPUTS="${IND_INPUTS:-${STRICT_ROOT}/16_point_shared_independent_repartition/01_inputs}"
SHARED_FIT="${SHARED_FIT:-${STRICT_ROOT}/25_point_fixed_ab_physical_tracks/01_model_run/normalized/point_model_fit.tsv}"
OUT_ROOT="${OUT_ROOT:-${STRICT_ROOT}/26_point_partition_models_fixedAB}"
POINT_BIN="${POINT_BIN:-${ROOT_DIR}/software/POInT/POInT_fixedAB}"
MODEL_FILE="${MODEL_FILE:-${ROOT_DIR}/software/POInT/models/WGD_bias_fix_fixedAB_model.txt}"
SUMMARIZE="${SUMMARIZE:-${ROOT_DIR}/scripts/fig2_summarize_point_partition_models.py}"
PLOT="${PLOT:-${ROOT_DIR}/scripts/plot_fig2_panel_d_point_partition_models.py}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
POINT_TOL="${POINT_TOL:-1e-7}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${SHARED_FIT}"
require_file "${MODEL_FILE}"
require_file "${SUMMARIZE}"
require_file "${PLOT}"
[[ -x "${POINT_BIN}" ]] || die "POInT_fixedAB binary unavailable: ${POINT_BIN}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT"

PAIR_ROOT="${OUT_ROOT}/01_pair_inputs"
IND_ROOT="${OUT_ROOT}/01_independent_inputs"
RUN_ROOT="${OUT_ROOT}/02_pair_runs"
IND_RUN_ROOT="${OUT_ROOT}/02_independent_runs"
SUMMARY_ROOT="${OUT_ROOT}/03_summary"
PLOT_ROOT="${OUT_ROOT}/04_panel_d"
mkdir -p "${RUN_ROOT}" "${IND_RUN_ROOT}" "${SUMMARY_ROOT}" "${PLOT_ROOT}"

# ---- copy inputs (POInT writes the .out tree next to the tree file) ----
cp -r "${PAIR_INPUTS}" "${PAIR_ROOT}"
cp -r "${IND_INPUTS}" "${IND_ROOT}"
rm -f "${PAIR_ROOT}"/*/species_tree.nex.out "${IND_ROOT}"/*/species_tree.nex.out

N_PILLARS=$(awk 'NR > 1 && NF {n++} END {print n + 0}' "${PAIR_ROOT}/PZI_PSP/pillars.txt")
echo "[INFO] Pillars: ${N_PILLARS}"

# ---- three pairwise shared fits (fixedAB) ----
for pair in PZI_PSP PZI_PPA PSP_PPA; do
    pdir="${PAIR_ROOT}/${pair}"
    rdir="${RUN_ROOT}/${pair}"
    mkdir -p "${rdir}"
    left="${pair%_*}"
    right="${pair#*_}"
    echo "[INFO] Pair fit (fixedAB): ${pair}"
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
            -fixedAB \
            > "${rdir}/point.stdout" \
            2> "${rdir}/point.stderr"
    )
    require_file "${pdir}/species_tree.nex.out"
    cp "${pdir}/species_tree.nex.out" "${rdir}/species_tree.nex.out"
done

# ---- three species-specific independent fits (fixedAB) ----
for species in PZI PSP PPA; do
    sdir="${IND_ROOT}/${species}"
    rdir="${IND_RUN_ROOT}/${species}"
    mkdir -p "${rdir}"
    echo "[INFO] Independent fit (fixedAB): ${species}"
    (
        cd "${rdir}"
        "${POINT_BIN}" \
            -d:2 \
            -g:"${sdir}/${species}_R1_geneorder.txt" \
            -g:"${sdir}/${species}_R2_geneorder.txt" \
            -g:"${sdir}/${species}_R3_geneorder.txt" \
            -o:"${sdir}/pillars.txt" \
            -m:"${MODEL_FILE}" \
            -t:"${sdir}/species_tree.nex" \
            -p:"${rdir}/posterior_probs.txt" \
            -z -l:"${POINT_TOL}" \
            -fixedAB \
            > "${rdir}/point.stdout" \
            2> "${rdir}/point.stderr"
    )
    require_file "${sdir}/species_tree.nex.out"
    cp "${sdir}/species_tree.nex.out" "${rdir}/species_tree.nex.out"
done

# ---- five-model summary and panel ----
"${PYTHON_BIN}" "${SUMMARIZE}" \
    --shared-fit "${SHARED_FIT}" \
    --pair-root "${RUN_ROOT}" \
    --independent-root "${IND_RUN_ROOT}" \
    --pillars "${N_PILLARS}" \
    --out "${SUMMARY_ROOT}/five_model_comparison.tsv"

"${PYTHON_BIN}" "${PLOT}" \
    --summary "${SUMMARY_ROOT}/five_model_comparison.tsv" \
    --outdir "${PLOT_ROOT}"

cat > "${OUT_ROOT}/software_and_parameters.txt" <<EOF
method=modified POInT v1.63 with fixed subgenome identity (-fixedAB build)
model=WGD_bias_fix_fixedAB_model.txt
pillars=${N_PILLARS}
shared_all=rooted three-taxon shared fit with physical A/B columns fixed
pair_shared=two-taxon rooted shared fit with physical A/B columns fixed plus independent singleton likelihood
all_independent=three species-specific replicated-terminal fits with physical A/B columns fixed; effective lnL is total/3
pair_parameter_count=3 process + 3 rooted two-taxon branches + 4 independent singleton = 10
all_shared_parameter_count=3 process + 5 rooted three-taxon branches = 8
all_independent_parameter_count=3 species x (3 process + 1 effective terminal branch) = 12
bic_sample_unit=common ${N_PILLARS} pillars for all five models
EOF
printf 'FixedAB five-model POInT partition analysis completed.\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] FixedAB five-model comparison completed: ${OUT_ROOT}"
