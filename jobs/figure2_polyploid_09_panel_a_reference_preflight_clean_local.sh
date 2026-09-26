#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Run the strict Figure 2a preflight with the reproducible WGD v2 environment.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
ENV_NAME="${ENV_NAME:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"

die() { echo "[ERROR] $*" >&2; exit 2; }

ENV_BIN="${CONDA_BASE}/envs/${ENV_NAME}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
[[ -x "${ENV_BIN}/python3" ]] || die "Missing clean WGD Python: ${ENV_BIN}/python3"
[[ -x "${TOOLS_BIN}/blastp" ]] || die "Missing tools environment executable: ${TOOLS_BIN}/blastp"

export PATH="${ENV_BIN}:${TOOLS_BIN}:${PATH}"
exec bash "${ROOT_DIR}/jobs/figure2_polyploid_09_panel_a_reference_preflight_local.sh"
