#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/reference_anchor_states}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"

INPUT_DIR="${KS_ROOT}/jcvi/inputs"
PAIR_SCRIPT="${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py"
STATE_SCRIPT="${ROOT_DIR}/scripts/fig2_polyploid_anchor_states.py"
PAIR_DIR="${OUT_ROOT}/pairs"
SUMMARY_DIR="${OUT_ROOT}/summary"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

[[ -e "${KS_ROOT}/jcvi/COMPLETE" ]] || die "Missing completion marker: ${KS_ROOT}/jcvi/COMPLETE"
[[ ! -e "${OUT_ROOT}/COMPLETE" ]] || die "Output is already complete: ${OUT_ROOT}" 
require_file "${PAIR_SCRIPT}"
require_file "${STATE_SCRIPT}"
for label in pzi_hapA1 pzi_hapB1 pzi_hapA2 pzi_hapB2 psp_hapA1 psp_hapB1 psp_hapA2 psp_hapB2 gtr; do
    for ext in bed pep; do
        require_file "${INPUT_DIR}/${label}/${label}.${ext}"
    done
done

mkdir -p "${OUT_ROOT}/logs" "${PAIR_DIR}" "${SUMMARY_DIR}"

set +u
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
set -u
command -v diamond >/dev/null 2>&1 || die "diamond unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi unavailable in ${JCVI_ENV}"

cat > "${SUMMARY_DIR}/reference_pair_plan.tsv" <<'EOF'
comparison	query_label	ref_label	comparison_type	analysis_role
pzi_hapA1_vs_gtr	pzi_hapA1	gtr	gtr_reference	main
pzi_hapB1_vs_gtr	pzi_hapB1	gtr	gtr_reference	main
psp_hapA1_vs_gtr	psp_hapA1	gtr	gtr_reference	main
psp_hapB1_vs_gtr	psp_hapB1	gtr	gtr_reference	main
pzi_hapA2_vs_gtr	pzi_hapA2	gtr	gtr_reference	sensitivity
pzi_hapB2_vs_gtr	pzi_hapB2	gtr	gtr_reference	sensitivity
psp_hapA2_vs_gtr	psp_hapA2	gtr	gtr_reference	sensitivity
psp_hapB2_vs_gtr	psp_hapB2	gtr	gtr_reference	sensitivity
EOF

while IFS=$'\t' read -r comparison query ref comparison_type analysis_role; do
    [[ "${comparison}" == "comparison" ]] && continue
    work="${PAIR_DIR}/${comparison}"
    mkdir -p "${work}"
    for label in "${query}" "${ref}"; do
        for ext in bed pep; do
            ln -sfn "${INPUT_DIR}/${label}/${label}.${ext}" "${work}/${label}.${ext}"
        done
    done
    cd "${work}"
    raw="${query}.${ref}.anchors"
    screened="${query}.${ref}.anchors.new"
    if [[ ! -s "${raw}" ]]; then
        python -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
            --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" \
            --min_size="${MIN_ANCHORS}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    fi
    require_file "${work}/${raw}"
    if [[ ! -s "${screened}" ]]; then
        python -m jcvi.compara.synteny screen --qbed="${query}.bed" --sbed="${ref}.bed" \
            --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
    fi
    require_file "${work}/${screened}"
    cd "${ROOT_DIR}"
done < "${SUMMARY_DIR}/reference_pair_plan.tsv"

python "${PAIR_SCRIPT}" \
    --pair-plan "${SUMMARY_DIR}/reference_pair_plan.tsv" \
    --pair-root "${PAIR_DIR}" \
    --out-pairs "${SUMMARY_DIR}/gtr_syntenic_pairs.tsv" \
    --out-blocks "${SUMMARY_DIR}/gtr_syntenic_blocks.tsv" \
    --out-status "${SUMMARY_DIR}/gtr_pair_status.tsv"

python "${STATE_SCRIPT}" \
    --pairs "${SUMMARY_DIR}/gtr_syntenic_pairs.tsv" \
    --reference-label gtr \
    --group main,PZI,pzi_hapA1,pzi_hapB1 main,PSP,psp_hapA1,psp_hapB1 \
            sensitivity,PZI,pzi_hapA2,pzi_hapB2 sensitivity,PSP,psp_hapA2,psp_hapB2 \
    --out-states "${SUMMARY_DIR}/gtr_anchor_states.tsv" \
    --out-summary "${SUMMARY_DIR}/gtr_anchor_state_summary.tsv" \
    --out-joint "${SUMMARY_DIR}/gtr_joint_anchor_states.tsv" \
    --out-candidates "${SUMMARY_DIR}/gtr_shared_reduced_anchor_candidates.tsv"

for output in \
    "${SUMMARY_DIR}/gtr_syntenic_pairs.tsv" \
    "${SUMMARY_DIR}/gtr_pair_status.tsv" \
    "${SUMMARY_DIR}/gtr_anchor_states.tsv" \
    "${SUMMARY_DIR}/gtr_anchor_state_summary.tsv" \
    "${SUMMARY_DIR}/gtr_joint_anchor_states.tsv" \
    "${SUMMARY_DIR}/gtr_shared_reduced_anchor_candidates.tsv"; do
    require_file "${output}"
done

{
    echo "python=$(python --version 2>&1)"
    echo "jcvi_env=${JCVI_ENV}"
    echo "threads=${THREADS}"
    echo "cscore=${CSCORE}"
    echo "min_anchors=${MIN_ANCHORS}"
    echo "minspan=${MINSPAN}"
    echo "reference=gtr"
    echo "interpretation=unanchored_loci_are_not_gene_losses"
} > "${OUT_ROOT}/software_and_parameters.txt"

touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Reference-anchor states completed: ${OUT_ROOT}"
