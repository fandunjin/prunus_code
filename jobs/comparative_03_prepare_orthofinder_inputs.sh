#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/comparative_genomics/orthofinder_input" "${ROOT_DIR}/logs"

python3 "${ROOT_DIR}/scripts/prepare_orthofinder_pep_inputs.py" "${ROOT_DIR}"

echo "OrthoFinder protein inputs:"
find "${ROOT_DIR}/comparative_genomics/orthofinder_input/pep" -maxdepth 1 -type f -name "*.pep.fa" | sort
echo
cat "${ROOT_DIR}/comparative_genomics/orthofinder_input/pep_manifest.tsv"
