#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SOURCE_DIR="${SOURCE_DIR:-${ROOT_DIR}/0903/haploid_repartition/inputs}"
BASE_DIR="${BASE_DIR:-${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs/main/pep}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs_copyAB14/pep}"
MANIFEST="${MANIFEST:-${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs_copyAB14/pep_manifest.tsv}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_copyab14_inputs.py" \
  --source "${SOURCE_DIR}" \
  --base "${BASE_DIR}" \
  --out "${OUT_DIR}" \
  --manifest "${MANIFEST}"
