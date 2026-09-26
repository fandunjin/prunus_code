#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Validate GTR-centered lower pair families using WGD v2's exact CDS rules.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"

die() { echo "[ERROR] $*" >&2; exit 2; }

PYTHON_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
PREP="${BASE}/01_prepared_inputs/lower"
TRACKS="${BASE}/../04_fractbias_8track_reference_exact/01_catalog/tracks"
OUT_ROOT="${BASE}/06_cds_preflight_lower"
VALIDATOR="${ROOT_DIR}/scripts/fig2_validate_wgd_pair_families.py"

[[ -s "${VALIDATOR}" ]] || die "Missing validator: ${VALIDATOR}"
mkdir -p "${OUT_ROOT}"

for track in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B PCAM PPE; do
  families="${PREP}/${track}.families.tsv"
  query_cds="${TRACKS}/${track}.cds"
  gtr_cds="${TRACKS}/GTR.cds"
  outdir="${OUT_ROOT}/${track}"
  [[ -s "${families}" ]] || die "Missing lower family table: ${families}"
  [[ -s "${query_cds}" ]] || die "Missing query CDS: ${query_cds}"
  [[ -s "${gtr_cds}" ]] || die "Missing GTR CDS: ${gtr_cds}"
  if [[ -e "${outdir}/COMPLETE" ]]; then
    echo "[OK] Reusing completed CDS validation: ${outdir}"
    continue
  fi
  "${PYTHON_BIN}" "${VALIDATOR}" \
    --families "${families}" \
    --cds "${gtr_cds}" \
    --cds "${query_cds}" \
    --dataset "${track}_GTR" \
    --outdir "${outdir}"
done

{
  head -n 1 "${OUT_ROOT}/PZI_A/PZI_A_GTR.cds_validation_summary.tsv"
  for track in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B PCAM PPE; do
    tail -n 1 "${OUT_ROOT}/${track}/${track}_GTR.cds_validation_summary.tsv"
  done
} > "${OUT_ROOT}/lower_cds_validation_summary.tsv"

printf 'Lower GTR-centered WGD CDS validation completed.\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] Figure 2a lower CDS validation completed: ${OUT_ROOT}"
cat "${OUT_ROOT}/lower_cds_validation_summary.tsv"
