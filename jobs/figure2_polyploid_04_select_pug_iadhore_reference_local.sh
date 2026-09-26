#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Reference-adapted Fig. 2b family screen using i-ADHoRe A/B pairs.

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
FROZEN="${FROZEN:-${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main}"
ORTHOGROUPS="${ORTHOGROUPS:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13/Orthogroups/Orthogroups.tsv}"
IADHORE="${IADHORE:-${BASE}/00_iadhore_reference_exact/02_ab_pairs}"
OUT="${OUT:-${BASE}/00_family_selection_reference_iadhore_allastral}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_file "${ROOT_DIR}/scripts/fig2_select_wgd_collinear_families.py"
require_file "${ORTHOGROUPS}"
[[ ! -e "${OUT}" ]] || die "Output directory already exists: ${OUT}"
require_file "${IADHORE}/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors.new"
require_file "${IADHORE}/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors.new"
require_file "${IADHORE}/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors.new"

ALL_LABELS=(pzi_A pzi_B psp_A psp_B ppa_A ppa_B pmal pav pspe pcam psal par pmu pdu ppe gtr rch)
SINGLE_COPY=(pmal pav pspe pcam psal par pmu pdu ppe)
for label in "${ALL_LABELS[@]}"; do
    require_file "${FROZEN}/id_maps/${label}.id_map.tsv"
    require_file "${FROZEN}/pep/${label}.fa"
done

cmd=(python3 "${ROOT_DIR}/scripts/fig2_select_wgd_collinear_families.py"
    --orthogroups "${ORTHOGROUPS}" --outdir "${OUT}"
    --outgroup gtr --second-outgroup rch --second-outgroup-mode at_least_one
    --collinear-pair "pzi_A,pzi_B,${IADHORE}/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors.new"
    --collinear-pair "psp_A,psp_B,${IADHORE}/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors.new"
    --collinear-pair "ppa_A,ppa_B,${IADHORE}/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors.new")
for label in "${ALL_LABELS[@]}"; do
    cmd+=(--id-map "${label}=${FROZEN}/id_maps/${label}.id_map.tsv")
    cmd+=(--pep "${label}=${FROZEN}/pep/${label}.fa")
done
for label in "${SINGLE_COPY[@]}"; do cmd+=(--required-single-copy "${label}"); done
"${cmd[@]}"

{
    echo "selection=reference_adapted_all_ASTRAL_taxa"
    echo "collinearity=i-ADHoRe_3.0_anchorpoints"
    echo "duplicate_evidence=PZI_PSP_PPA_unique_AB_collinear_pairs"
    echo "single_copy_background=${SINGLE_COPY[*]}"
    echo "removed_constraints=PCAM_GTR_collinearity;PPE_single_copy_only"
    echo "JCVI_anchor_inputs=forbidden"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] i-ADHoRe reference family selection completed: ${OUT}"
