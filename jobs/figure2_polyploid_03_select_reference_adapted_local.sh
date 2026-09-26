#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Reference-adapted Fig. 2b family screen with all taxa in the ASTRAL tree.
# Only PZI/PSP/PPA A-B collinear pairs are used as duplicate evidence.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
FROZEN="${FROZEN:-${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main}"
ORTHOGROUPS="${ORTHOGROUPS:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13/Orthogroups/Orthogroups.tsv}"
PAIRS="${PAIRS:-${ROOT_DIR}/downstream_analysis/ks/jcvi/pairs}"
OUT="${OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/00_family_selection_reference_adapted_allastral}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${ROOT_DIR}/scripts/fig2_select_wgd_collinear_families.py"
require_file "${ORTHOGROUPS}"
[[ ! -e "${OUT}" ]] || die "Output directory already exists: ${OUT}"
command -v python3 >/dev/null 2>&1 || die "python3 is unavailable on PATH"

# These labels are the terminal taxa represented by the current ASTRAL tree.
ALL_LABELS=(pzi_A pzi_B psp_A psp_B ppa_A ppa_B pmal pav pspe pcam psal par pmu pdu ppe gtr rch)
SINGLE_COPY=(pmal pav pspe pcam psal par pmu pdu ppe)

for label in "${ALL_LABELS[@]}"; do
    require_file "${FROZEN}/id_maps/${label}.id_map.tsv"
    require_file "${FROZEN}/pep/${label}.fa"
done

require_file "${PAIRS}/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors.new"
require_file "${PAIRS}/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors.new"
require_file "${PAIRS}/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors.new"

cmd=(python3 "${ROOT_DIR}/scripts/fig2_select_wgd_collinear_families.py"
    --orthogroups "${ORTHOGROUPS}"
    --outdir "${OUT}"
    --outgroup gtr
    --second-outgroup rch
    --second-outgroup-mode at_least_one
    --collinear-pair "pzi_A,pzi_B,${PAIRS}/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors.new"
    --collinear-pair "psp_A,psp_B,${PAIRS}/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors.new"
    --collinear-pair "ppa_A,ppa_B,${PAIRS}/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors.new")
for label in "${ALL_LABELS[@]}"; do
    cmd+=(--id-map "${label}=${FROZEN}/id_maps/${label}.id_map.tsv")
    cmd+=(--pep "${label}=${FROZEN}/pep/${label}.fa")
done
for label in "${SINGLE_COPY[@]}"; do
    cmd+=(--required-single-copy "${label}")
done
"${cmd[@]}"

require_file "${OUT}/selected_collinear_wgd_families.tsv"
require_file "${OUT}/screening_summary.tsv"
[[ -d "${OUT}/family_fastas" ]] || die "No selected family FASTAs were written"
{
    echo "selection=reference_adapted_all_ASTRAL_taxa"
    echo "duplicate_evidence=PZI_PSP_PPA_unique_AB_collinear_pairs"
    echo "single_copy_background=${SINGLE_COPY[*]}"
    echo "outgroup=gtr_exactly_one"
    echo "second_outgroup=rch_at_least_one"
    echo "removed_constraints=PCAM_GTR_collinearity;PPE_single_copy_only;RCH_exactly_one"
    echo "orthogroups=${ORTHOGROUPS}"
    echo "analysis_scope=primary_A1_B1_only"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] Reference-adapted all-ASTRAL-taxa family selection completed: ${OUT}"
