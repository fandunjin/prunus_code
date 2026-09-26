#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
ALIGN="${BASE}/04_alignments"
OUT="${OUT:-${BASE}/05_gene_trees}"
IQTREE_BIN="${IQTREE_BIN:-}"
THREADS="${THREADS:-8}"
UFB_REPLICATES="${UFB_REPLICATES:-2000}"
ALRT_REPLICATES="${ALRT_REPLICATES:-1000}"
SEED="${SEED:-20260825}"

[[ -e "${ALIGN}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_03_align_trim_local.sh first" >&2; exit 2; }
mkdir -p "${OUT}"
if [[ -z "${IQTREE_BIN}" ]]; then
  for candidate in iqtree3 iqtree2 iqtree; do command -v "${candidate}" >/dev/null 2>&1 && IQTREE_BIN="$(command -v "${candidate}")" && break; done
fi
[[ -n "${IQTREE_BIN}" ]] || { echo "[ERROR] IQ-TREE unavailable; set IQTREE_BIN" >&2; exit 2; }
run_one() {
  input="$1"; dataset="$2"; base=$(basename "${input}" .fa)
  if [[ -s "${OUT}/${dataset}/genes/${base}.treefile" && -s "${OUT}/${dataset}/genes/${base}.ufboot" && $(wc -l < "${OUT}/${dataset}/genes/${base}.ufboot") -ge "${UFB_REPLICATES}" ]]; then
    return 0
  fi
  "${IQTREE_BIN}" -s "${input}" -m MFP -B "${UFB_REPLICATES}" --alrt "${ALRT_REPLICATES}" --bnni -T 1 -seed "${SEED}" --prefix "${OUT}/${dataset}/genes/${base}" > "${OUT}/${dataset}/logs/${base}.iqtree.log" 2>&1
}
export IQTREE_BIN UFB_REPLICATES ALRT_REPLICATES SEED OUT
export -f run_one
for dataset in main sensitivity; do
  mkdir -p "${OUT}/${dataset}/genes" "${OUT}/${dataset}/logs"
  expected=$(wc -l < "${ALIGN}/${dataset}/qc/retained_orthogroups.txt")
  find "${ALIGN}/${dataset}/qc/retained" -maxdepth 1 -type f -name '*.fa' -print0 | sort -z | \
    xargs -0 -r -n 1 -P "${THREADS}" bash -c 'run_one "$1" "'"${dataset}"'"' _
  observed=$(find "${OUT}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -size +0c | wc -l)
  bootstrapped=$(find "${OUT}/${dataset}/genes" -maxdepth 1 -type f -name '*.ufboot' -size +0c | wc -l)
  [[ "${observed}" -eq "${expected}" && "${bootstrapped}" -eq "${expected}" ]] || { echo "[ERROR] ${dataset}: expected ${expected}, observed ${observed}/${bootstrapped}" >&2; exit 2; }
done
"${IQTREE_BIN}" --version > "${OUT}/software_versions.txt" 2>&1
printf 'seed=%s\nufboot_replicates=%s\nsh_alrt_replicates=%s\n' "${SEED}" "${UFB_REPLICATES}" "${ALRT_REPLICATES}" >> "${OUT}/software_versions.txt"
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned gene trees completed: ${OUT}"
