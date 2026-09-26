#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Validate upper A-B pair families using WGD v2's exact Biopython CDS rules.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"

die() { echo "[ERROR] $*" >&2; exit 2; }

PYTHON_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
MCL_ROOT="${BASE}/03_blastp_mcl"
CDS_ROOT="${BASE}/01_prepared_inputs/upper"
OUT_ROOT="${BASE}/04_cds_preflight_upper"
VALIDATOR="${ROOT_DIR}/scripts/fig2_validate_wgd_pair_families.py"

[[ -s "${VALIDATOR}" ]] || die "Missing validator: ${VALIDATOR}"
mkdir -p "${OUT_ROOT}"

for dataset in PZI PSP PPA; do
  families="${MCL_ROOT}/${dataset}/filtered_pairs/${dataset}_AB.families.tsv"
  cds="${CDS_ROOT}/${dataset}_AB.cds.fa"
  outdir="${OUT_ROOT}/${dataset}"
  [[ -s "${families}" ]] || die "Missing MCL-filtered families: ${families}"
  [[ -s "${cds}" ]] || die "Missing CDS FASTA: ${cds}"
  if [[ -e "${outdir}/COMPLETE" ]]; then
    echo "[OK] Reusing completed CDS validation: ${outdir}"
    continue
  fi
  "${PYTHON_BIN}" "${VALIDATOR}" \
    --families "${families}" \
    --cds "${cds}" \
    --dataset "${dataset}_AB" \
    --outdir "${outdir}"
done

{
  head -n 1 "${OUT_ROOT}/PZI/PZI_AB.cds_validation_summary.tsv"
  for dataset in PZI PSP PPA; do
    tail -n 1 "${OUT_ROOT}/${dataset}/${dataset}_AB.cds_validation_summary.tsv"
  done
} > "${OUT_ROOT}/upper_cds_validation_summary.tsv"

printf 'Upper A-B WGD CDS validation completed.\n' > "${OUT_ROOT}/COMPLETE"
echo "[OK] Figure 2a upper CDS validation completed: ${OUT_ROOT}"
cat "${OUT_ROOT}/upper_cds_validation_summary.tsv"
