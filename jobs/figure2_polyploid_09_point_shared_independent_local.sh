#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Compare a common ancestral WGD with three independent WGD fits.
# The independent fits use three identical POInT terminals per species because
# POInT 1.63 requires a rooted multi-terminal input. Root/internal branches are
# fixed at zero; the replicated lnL is divided by three before AIC/BIC.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SHARED_ROOT="${SHARED_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_official_reference_exact_r13}"
PILLARS="${PILLARS:-${SHARED_ROOT}/02_ancestral_order/point_pillars_optimized.txt}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_point_shared_vs_independent}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
POINT_BIN="${POINT_BIN:-${ROOT_DIR}/software/POInT/POInT}"
MODEL_FILE="${MODEL_FILE:-${ROOT_DIR}/software/POInT/models/WGD_bias_fix_model.txt}"
PREPARE="${ROOT_DIR}/scripts/fig2_prepare_point_independent_replicates.py"
SUMMARIZE="${ROOT_DIR}/scripts/fig2_summarize_point_shared_independent.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${SHARED_ROOT}/04_summary/point_model_fit_all.tsv"
require_file "${PILLARS}"
require_file "${MODEL_FILE}"
require_file "${PREPARE}"
require_file "${SUMMARIZE}"
[[ -x "${POINT_BIN}" ]] || die "POInT binary unavailable: ${POINT_BIN}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT"

mkdir -p "${OUT_ROOT}"
INPUT_ROOT="${OUT_ROOT}/01_inputs"
RUN_ROOT="${OUT_ROOT}/02_independent_runs"
mkdir -p "${RUN_ROOT}"

"${PYTHON_BIN}" "${PREPARE}" \
    --inputs "${SHARED_ROOT}/01_inputs" \
    --optimized "${PILLARS}" \
    --outdir "${INPUT_ROOT}"

for species in PZI PSP PPA; do
    sdir="${INPUT_ROOT}/${species}"
    rdir="${RUN_ROOT}/${species}"
    mkdir -p "${rdir}"
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
            -z -l:1e-7 \
            > "${rdir}/point.stdout" \
            2> "${rdir}/point.stderr"
    )
    # POInT writes the .out tree next to the tree file, not next to stdout.
    require_file "${sdir}/species_tree.nex.out"
    cp "${sdir}/species_tree.nex.out" "${rdir}/species_tree.nex.out"
done

${PYTHON_BIN} "${SUMMARIZE}" \
    --shared-fit "${SHARED_ROOT}/04_summary/point_model_fit_all.tsv" \
    --independent-root "${RUN_ROOT}" \
    --pillars "$(($(wc -l < "${PILLARS}") - 1))" \
    --out "${OUT_ROOT}/03_shared_vs_independent.tsv"

cat > "${OUT_ROOT}/software_and_parameters.txt" <<EOF
method=official POInT v1.63
shared_model=existing reference_WGD_bias_fix fit on (PPA,(PSP,PZI))
independent_model=three species-specific replicated-terminal fits; effective lnL is total/3
species=PZI,PSP,PPA
pillars=$(($(wc -l < "${PILLARS}") - 1))
model=WGD_bias_fix_model.txt (3 process parameters)
independent_parameter_count=3 species x (3 process parameters + 1 effective terminal branch) = 12
shared_parameter_count=3 process parameters + 5 POInT branches (including shared root/stem) = 8
root_state=Dupl
zero_branch_option=-z
optimization_tolerance=1e-7
interpretation=This is a likelihood comparison of a shared ancestral WGD process versus separate species-specific WGD processes; it is not evidence of parental identity.
EOF
printf 'Shared versus independent POInT comparison completed.\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] Shared/independent POInT analysis completed: ${OUT_ROOT}"
