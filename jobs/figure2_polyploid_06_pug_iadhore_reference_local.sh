#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Reference-method PUG mapping for the i-ADHoRe-derived family trees.

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_reference_iadhore_allastral}"
TREE_STAGE="${TREE_STAGE:-${BASE}/02_family_trees_reference_iadhore_allastral}"
OUT="${OUT:-${BASE}/03_pug_reference_iadhore_allastral}"
SPECIES_TREE_FILE="${SPECIES_TREE_FILE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/reference_adapted_pug_species_tree.nwk}"
[[ -s "${SPECIES_TREE_FILE}" ]] || { echo "[ERROR] Missing species tree: ${SPECIES_TREE_FILE}" >&2; exit 2; }
PUG_SPECIES_TREE="$(tr -d '\r\n' < "${SPECIES_TREE_FILE}")"
PUG_EXPECTED_LABELS="pzi_A,pzi_B,psp_A,psp_B,ppa_A,ppa_B,pmal,pav,pspe,pcam,psal,par,pmu,pdu,ppe,gtr,rch"
export ROOT_DIR INPUT TREE_STAGE OUT PUG_SPECIES_TREE PUG_EXPECTED_LABELS
exec bash "${ROOT_DIR}/jobs/figure2_polyploid_05_pug_local.sh"
