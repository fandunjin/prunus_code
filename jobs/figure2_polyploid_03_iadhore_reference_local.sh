#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Reference-method Fig. 2b collinearity: i-ADHoRe 3.0 on orthogroup families.
# The resulting A/B anchor pairs feed the PUG family screen; JCVI anchors are
# intentionally not accepted here.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CATALOG_DIR="${CATALOG_DIR:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/01_catalog}"
ORTHOGROUPS="${ORTHOGROUPS:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13/Orthogroups/Orthogroups.tsv}"
ID_MAPS="${ID_MAPS:-${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/id_maps}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/00_iadhore_reference_exact}"
THREADS="${THREADS:-12}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IADHORE_BIN="${IADHORE_BIN:-${ROOT_DIR}/software/i-ADHoRe/install/bin/i-adhore}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_exec() { [[ -x "$1" ]] || die "Required executable unavailable: $1"; }

[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}"
require_file "${ROOT_DIR}/scripts/fig2_prepare_iadhore_orthofinder.py"
require_file "${ROOT_DIR}/scripts/fig2_extract_iadhore_ab_pairs.py"
require_file "${ORTHOGROUPS}"
for track in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B; do require_file "${CATALOG_DIR}/tracks/${track}.bed"; done
require_exec "${IADHORE_BIN}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"

PREP_DIR="${OUT_ROOT}/01_input"
mkdir -p "${OUT_ROOT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_prepare_iadhore_orthofinder.py" \
    --catalog "${CATALOG_DIR}" --orthogroups "${ORTHOGROUPS}" \
    --id-maps "${ID_MAPS}" --outdir "${PREP_DIR}" --threads "${THREADS}"
require_file "${PREP_DIR}/iadhore.ini"
mkdir -p "${PREP_DIR}/iadhore_output"
"${IADHORE_BIN}" "${PREP_DIR}/iadhore.ini" \
    > "${OUT_ROOT}/iadhore.stdout.log" 2> "${OUT_ROOT}/iadhore.stderr.log"
require_file "${PREP_DIR}/iadhore_output/multiplicons.txt"
require_file "${PREP_DIR}/iadhore_output/anchorpoints.txt"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_extract_iadhore_ab_pairs.py" \
    --iadhore-output "${PREP_DIR}/iadhore_output" \
    --stable-map "${PREP_DIR}/stable_gene_map.tsv" \
    --outdir "${OUT_ROOT}/02_ab_pairs"
require_file "${OUT_ROOT}/02_ab_pairs/COMPLETE"

{
    echo "method=i-ADHoRe 3.0 gene-order family collinearity"
    echo "family_table=${ORTHOGROUPS} translated to stable gene IDs"
    echo "gap_size=30"
    echo "cluster_gap=35"
    echo "q_value=0.75"
    echo "prob_cutoff=0.01"
    echo "anchor_points=3"
    echo "alignment_method=gg2"
    echo "level_2_only=false"
    echo "multiple_hypothesis_correction=FDR"
    echo "PZI_PSP_PPA_A_B_pairs=${OUT_ROOT}/02_ab_pairs"
    echo "JCVI_anchor_inputs=forbidden"
} > "${OUT_ROOT}/software_and_parameters.txt"
{
    "${IADHORE_BIN}" --version 2>&1 | head -1 || true
    "${PYTHON_BIN}" --version 2>&1 | head -1
} > "${OUT_ROOT}/tool_versions.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Reference i-ADHoRe collinearity completed: ${OUT_ROOT}"
