#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
SOURCE="${BASE}/03_single_copy"
INPUTS="${BASE}/01_frozen_inputs"
OUT="${OUT:-${BASE}/04_alignments}"
THREADS="${THREADS:-12}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[[ -e "${SOURCE}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_02_extract_single_copy_local.sh first" >&2; exit 2; }
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 2; }
for exe in mafft trimal; do command -v "${exe}" >/dev/null 2>&1 || { echo "[ERROR] ${exe} unavailable" >&2; exit 2; }; done
mkdir -p "${OUT}"

run_one() {
  input="$1"; dataset="$2"
  base=$(basename "${input}" .fa)
  mafft --localpair --maxiterate 1000 --thread 1 "${input}" > "${OUT}/${dataset}/mafft/${base}.fa" 2> "${OUT}/${dataset}/logs/${base}.mafft.log"
  trimal -in "${OUT}/${dataset}/mafft/${base}.fa" -out "${OUT}/${dataset}/trimmed/${base}.fa" -fasta -automated1 > "${OUT}/${dataset}/logs/${base}.trimal.log" 2>&1
}
export OUT
export -f run_one
for dataset in main sensitivity; do
  mkdir -p "${OUT}/${dataset}/mafft" "${OUT}/${dataset}/trimmed" "${OUT}/${dataset}/logs"
  find "${SOURCE}/${dataset}/unaligned" -maxdepth 1 -type f -name '*.fa' -print0 | sort -z | \
    xargs -0 -r -n 1 -P "${THREADS}" bash -c 'run_one "$1" "'"${dataset}"'"' _
  find "${INPUTS}/${dataset}/pep" -maxdepth 1 -type f -name '*.fa' -printf '%f\n' | sed 's/\.fa$//' | sort > "${OUT}/${dataset}/species.txt"
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_alignment_qc.py" \
    --alignment-dir "${OUT}/${dataset}/trimmed" --species-file "${OUT}/${dataset}/species.txt" \
    --out-dir "${OUT}/${dataset}/qc" --min-length 100 --min-informative 10 --max-missing 0.20 --min-alignments 200
done
{
  mafft --version
  trimal --version
} > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned alignments and QC completed: ${OUT}"
