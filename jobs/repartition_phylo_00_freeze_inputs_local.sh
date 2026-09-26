#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
MANIFEST="${MANIFEST:-${BASE}/config/input_manifest.tsv}"
OUT="${OUT:-${BASE}/01_frozen_inputs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SOURCE_DIRECTORY="${SOURCE_DIRECTORY:-${ROOT_DIR}/haploid_repartition}"

[[ -s "${MANIFEST}" ]] || { echo "[ERROR] Missing manifest: ${MANIFEST}" >&2; exit 2; }
if [[ -e "${OUT}" ]]; then
  if [[ -e "${OUT}/COMPLETE" ]]; then
    echo "[ERROR] Output exists and is complete: ${OUT}; freeze is immutable" >&2
    exit 2
  fi
  backup="${OUT}.incomplete.$(date +%Y%m%d_%H%M%S)"
  mv -- "${OUT}" "${backup}"
  echo "[WARN] Incomplete freeze moved to ${backup}" >&2
fi
mkdir -p "${BASE}/logs"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_prepare_inputs.py" \
  --manifest "${MANIFEST}" \
  --out-dir "${OUT}" \
  --skip-invalid-proteins
{
  date -Is
  "${PYTHON_BIN}" --version
  sha256sum "${MANIFEST}" "${BASE}/config/focal_quartets.tsv"
  printf 'source_directory=%s\n' "${SOURCE_DIRECTORY}"
  excluded_count=$(awk 'NR > 1 {n++} END {print n+0}' "${OUT}/excluded_invalid_proteins.tsv")
  printf 'excluded_invalid_proteins=%s\n' "${excluded_count}"
} > "${OUT}/provenance.txt"
echo "[OK] New repartitioned phylogenetic inputs frozen: ${OUT}"
