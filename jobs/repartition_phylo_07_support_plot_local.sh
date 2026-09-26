#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
OUT="${OUT:-${BASE}/08_final_summary}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[[ -e "${BASE}/07_topology_frequencies/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_06_topology_local.sh first" >&2; exit 2; }
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 2; }
# The plotting program creates OUT itself and refuses a pre-existing directory.
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_merge_supports_plot.py" \
  --base "${BASE}" --manifest "${BASE}/config/input_manifest.tsv" --out-dir "${OUT}" \
  --datasets main sensitivity --outgroup rch --width 7.2 --height 7.0
for dataset in main sensitivity; do
  [[ -s "${OUT}/${dataset}/final_astral_support_tree.pdf" ]] || { echo "[ERROR] Missing support tree for ${dataset}" >&2; exit 2; }
done
printf 'outgroup=rch\nnode_label_order=bootstrap/PP/gCF/likelihood-sCF\n' > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned support trees and plots completed: ${OUT}"
