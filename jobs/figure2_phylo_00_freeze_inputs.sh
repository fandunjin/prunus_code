#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
MANIFEST="${BASE}/config/input_manifest.tsv"
OUT="${BASE}/01_frozen_inputs"
PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${BASE}/logs"
if [ -e "${OUT}" ]; then
    echo "[ERROR] Frozen input directory already exists: ${OUT}" >&2
    echo "[ERROR] Freeze is immutable; move the old directory aside before rerunning." >&2
    exit 3
fi

GTR="${ROOT}/data/gtr/Gtr.NCBI.prot.faa"
GTR_MD5="70ee906f4df33d739be312eabcaecebd"
if [ ! -s "${GTR}" ]; then
    echo "[ERROR] Missing Gillenia protein FASTA: ${GTR}" >&2
    exit 3
fi

observed_md5="$(md5sum "${GTR}" | awk '{print tolower($1)}')"
if [ "${observed_md5}" != "${GTR_MD5}" ]; then
    echo "[ERROR] Gillenia protein MD5 mismatch: ${observed_md5}" >&2
    echo "[ERROR] Expected ${GTR_MD5}; verify data/gtr/Gtr.NCBI.prot.faa before rerunning." >&2
    exit 3
fi

"${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_prepare_inputs.py" \
    --manifest "${MANIFEST}" \
    --out-dir "${OUT}"

{
    date -Is
    "${PYTHON_BIN}" --version
    sha256sum "${MANIFEST}" "${BASE}/config/focal_quartets.tsv"
} > "${OUT}/provenance.txt" 2>&1

echo "[OK] Immutable inputs: ${OUT}"
