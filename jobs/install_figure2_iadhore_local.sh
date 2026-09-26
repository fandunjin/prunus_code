#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Install only the tools required for Figure 2 self-synteny into an isolated Conda environment.

set -euo pipefail

ENV_NAME="${ENV_NAME:-figure2_iadhore}"
CONDA_EXE="${CONDA_EXE:-$(command -v conda || true)}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -n "${CONDA_EXE}" ]] || die "conda is unavailable on PATH"
CONDA_BASE="$(${CONDA_EXE} info --base)"
TARGET="${CONDA_BASE}/envs/${ENV_NAME}/bin"

if [[ -x "${TARGET}/i-adhore" && -x "${TARGET}/blastp" && -x "${TARGET}/makeblastdb" && -x "${TARGET}/mcl" ]]; then
  echo "[OK] Figure 2 i-ADHoRe tool environment already exists: ${TARGET}"
  exit 0
fi

"${CONDA_EXE}" create -y -n "${ENV_NAME}" --override-channels \
  -c conda-forge -c bioconda \
  'i-adhore=3.0.01' blast mcl

for tool in i-adhore blastp makeblastdb mcl; do
  [[ -x "${TARGET}/${tool}" ]] || die "Installation did not provide ${tool}: ${TARGET}/${tool}"
done
echo "[OK] Figure 2 i-ADHoRe tool environment: ${TARGET}"
