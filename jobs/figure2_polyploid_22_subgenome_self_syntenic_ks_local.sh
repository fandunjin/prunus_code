#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Independently calculate cross-chromosome self-syntenic Ks for six A/B tracks.
# This wrapper invokes the established one-dataset pipeline serially so each
# BLASTP/MCL/i-ADHoRe/CodeML analysis is insulated from the other five.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
THREADS="${THREADS:-12}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"

for dataset in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B; do
  echo "[INFO] Starting independent self-syntenic Ks: ${dataset}"
  ROOT_DIR="${ROOT_DIR}" THREADS="${THREADS}" WGD_ENV="${WGD_ENV}" TOOLS_ENV="${TOOLS_ENV}" \
    DATASET="${dataset}" bash "${ROOT_DIR}/jobs/figure2_polyploid_18_self_syntenic_ks_local.sh"
done

echo "[OK] All six independently run subgenome self-syntenic Ks analyses completed."
