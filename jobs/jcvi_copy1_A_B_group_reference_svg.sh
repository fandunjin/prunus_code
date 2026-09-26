#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

SCRIPT="${PROJECT_ROOT}/software/jcvi_karyotype_svg.pl"
ROOT="${HOME}/user/prunus/analyse/karyotype/copy1_anchor"
BED="${ROOT}/inputs"
ANCHOR="${ROOT}"
OUT="${OUT:-${ROOT}/svg_group_A_B_reference}"
LOG="${OUT}/logs"
SELF_BED="${HOME}/user/prunus/analyse/karyotype/copy1_anchor/inputs/pzi_copy1/pzi_copy1.bed"
SELF_ANCHORS="${HOME}/user/prunus/analyse/karyotype/copy1_anchor/pzi_copy1.pzi_copy1.anchors"
QUERY_BED="${BED}/pzi_copy1/pzi_copy1.bed"
MIN_BLOCK_SIZE="${MIN_BLOCK_SIZE:-1}"
MIN_ANCHOR_SCORE="${MIN_ANCHOR_SCORE:-0}"
MAX_QUERY_HITS="${MAX_QUERY_HITS:-0}"
MAX_SUBJECT_HITS="${MAX_SUBJECT_HITS:-0}"
UNIQUE_PAIRS="${UNIQUE_PAIRS:-0}"

mkdir -p "${OUT}" "${LOG}"

[[ -s "${SCRIPT}" ]] || { echo "Missing Perl renderer: ${SCRIPT}" >&2; exit 1; }
[[ -s "${QUERY_BED}" ]] || { echo "Missing PZI copy1 BED" >&2; exit 1; }
[[ -s "${BED}/pcam/pcam.bed" ]] || { echo "Missing PCAM BED" >&2; exit 1; }
[[ -s "${BED}/ppe/ppe.bed" ]] || { echo "Missing PPE BED" >&2; exit 1; }
[[ -s "${BED}/gtr/gtr.bed" ]] || { echo "Missing GTR BED" >&2; exit 1; }
[[ -s "${BED}/ppa/ppa.bed" ]] || { echo "Missing PPA BED" >&2; exit 1; }
[[ -s "${BED}/psp/psp.bed" ]] || { echo "Missing PSP BED" >&2; exit 1; }
[[ -s "${SELF_BED}" ]] || { echo "Missing copy1 self BED: ${SELF_BED}" >&2; exit 1; }
[[ -s "${SELF_ANCHORS}" ]] || { echo "Missing copy1 self anchors: ${SELF_ANCHORS}" >&2; exit 1; }

# Requested PZI copy1 order.  A+B below is the complete 16-chromosome copy1
# set, with the A group followed by the B group.
A="chr01,chr02,chr05,chr06,chr12,chr07,chr10,chr14"
B="chr13,chr15,chr11,chr16,chr03,chr04,chr08,chr09"
AB="chr01,chr02,chr05,chr06,chr12,chr07,chr10,chr14,chr13,chr15,chr11,chr16,chr03,chr04,chr08,chr09"

PCAM="CM114280.1,CM114281.1,CM114282.1,CM114283.1,CM114284.1,CM114285.1,CM114286.1,CM114287.1"
PPE="NC_034009.1,NC_034010.1,NC_034011.1,NC_034012.1,NC_034013.1,NC_034014.1,NC_034015.1,NC_034016.1"
GTR="Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09"
PPA="OZ060511.1,OZ060512.1,OZ060513.1,OZ060514.1,OZ060515.1,OZ060516.1,OZ060517.1,OZ060518.1,OZ060519.1,OZ060520.1,OZ060521.1,OZ060522.1,OZ060523.1,OZ060524.1,OZ060525.1,OZ060526.1"
PSP="chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12,chr13,chr14,chr15,chr16"

render() {
    local label="$1"
    local query_order="$2"
    local subject_order="$3"
    local subject_label="$4"
    local subject_bed="$5"
    local anchors="$6"
    local query_bed="$7"
    local group_name="$8"
    local colored_points="${9:-0}"
    local band_height="${10:-1}"
    local identity_diagonal="${11:-0}"
    local output="${OUT}/${label}.svg"
    local -a unique_args=()
    local -a point_args=()
    local -a identity_args=()

    if [[ "${UNIQUE_PAIRS}" == "1" ]]; then
        unique_args=(--unique-pairs)
    fi
    if [[ "${colored_points}" == "1" ]]; then
        point_args=(--colored-points)
    fi
    if [[ "${identity_diagonal}" == "1" ]]; then
        identity_args=(--identity-diagonal)
    fi

    [[ -s "${anchors}" ]] || { echo "Missing anchors: ${anchors}" >&2; exit 1; }
    echo "[RUN] ${label}"
    perl "${SCRIPT}" \
        --qbed "${query_bed}" \
        --sbed "${subject_bed}" \
        --anchors "${anchors}" \
        --query-seqids "${query_order}" \
        --subject-seqids "${subject_order}" \
        --mode bp \
        --min-block-size "${MIN_BLOCK_SIZE}" \
        --min-anchor-score "${MIN_ANCHOR_SCORE}" \
        --max-query-hits "${MAX_QUERY_HITS}" \
        --max-subject-hits "${MAX_SUBJECT_HITS}" \
        --band-height "${band_height}" \
        "${unique_args[@]}" \
        "${point_args[@]}" \
        "${identity_args[@]}" \
        --compact \
        --labels \
        --legend \
        --title "PZI copy1 ${group_name} vs ${subject_label}" \
        --query-label "PZI copy1 ${group_name} group" \
        --subject-label "${subject_label} chromosomes" \
        --output "${output}"
}

# A group versus PCAM, PPE and GTR.
render "pzi_copy1_A_vs_pcam" "${A}" "${PCAM}" "PCAM" \
    "${BED}/pcam/pcam.bed" "${ANCHOR}/pzi_copy1.pcam.lifted.anchors" "${QUERY_BED}" "A"
render "pzi_copy1_A_vs_ppe" "${A}" "${PPE}" "PPE" \
    "${BED}/ppe/ppe.bed" "${ANCHOR}/pzi_copy1.ppe.lifted.anchors" "${QUERY_BED}" "A"
render "pzi_copy1_A_vs_gtr" "${A}" "${GTR}" "GTR" \
    "${BED}/gtr/gtr.bed" "${ANCHOR}/pzi_copy1.gtr.lifted.anchors" "${QUERY_BED}" "A"

# B group versus PCAM, PPE and GTR.
render "pzi_copy1_B_vs_pcam" "${B}" "${PCAM}" "PCAM" \
    "${BED}/pcam/pcam.bed" "${ANCHOR}/pzi_copy1.pcam.lifted.anchors" "${QUERY_BED}" "B"
render "pzi_copy1_B_vs_ppe" "${B}" "${PPE}" "PPE" \
    "${BED}/ppe/ppe.bed" "${ANCHOR}/pzi_copy1.ppe.lifted.anchors" "${QUERY_BED}" "B"
render "pzi_copy1_B_vs_gtr" "${B}" "${GTR}" "GTR" \
    "${BED}/gtr/gtr.bed" "${ANCHOR}/pzi_copy1.gtr.lifted.anchors" "${QUERY_BED}" "B"

# Complete copy1 (A followed by B) versus PPA and PSP.
render "pzi_copy1_AB_vs_ppa" "${AB}" "${PPA}" "PPA" \
    "${BED}/ppa/ppa.bed" "${ANCHOR}/pzi_copy1.ppa.lifted.anchors" "${QUERY_BED}" "A+B"
render "pzi_copy1_AB_vs_psp" "${AB}" "${PSP}" "PSP" \
    "${BED}/psp/psp.bed" "${ANCHOR}/pzi_copy1.psp.lifted.anchors" "${QUERY_BED}" "A+B"

# Complete copy1 self-comparison, using the existing JCVI self anchors.
render "pzi_copy1_AB_self" "${AB}" "${AB}" "PZI copy1" \
    "${SELF_BED}" "${SELF_ANCHORS}" "${SELF_BED}" "A+B" 0 2 1

cat > "${OUT}/README.txt" <<EOF
PZI copy1 group-restricted JCVI SVG figures

A group:  ${A}
B group:  ${B}
A+B order: ${AB}

A/B groups were selected from the full pzi_copy1 BED and full lifted anchors;
the anchors were not recomputed.  PPE is restricted to NC_034009.1-NC_034016.1
(eight primary chromosomes), GTR to Chr01-Chr09, and PPA to OZ060511.1-OZ060526.1
(sixteen primary chromosomes).

All figures use base-pair mode. Filtering is controlled by MIN_BLOCK_SIZE,
MIN_ANCHOR_SCORE, MAX_QUERY_HITS, MAX_SUBJECT_HITS and UNIQUE_PAIRS.

The self-comparison uses the existing JCVI raw anchors from:
${SELF_ANCHORS}
EOF

echo "[OK] Outputs written to ${OUT}"
