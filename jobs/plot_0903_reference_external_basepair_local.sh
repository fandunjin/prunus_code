#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Render the 0903 seven-species JCVI ribbons locally in the reference-external
# style.  Run from the repository root in WSL or another bash environment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT="${ROOT}/scripts/plot_0903_reference_external_basepair.py"
OUT_DIR="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_multispecies_A1B1_pcam_pairwise"

[[ -s "${SCRIPT}" ]] || { echo "[ERROR] Missing renderer: ${SCRIPT}" >&2; exit 1; }
"${PYTHON_BIN}" "${SCRIPT}" \
  --hap-input "${ROOT}/0903/haploid_repartition/inputs" \
  --jcvi-input "${ROOT}/0903/downstream_analysis/figure2_repartition_ks/inputs" \
  --pair-work "${ROOT}/0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work" \
  --rch-input "${ROOT}/downstream_analysis/local_A1B1_ribbon_inputs" \
  --outdir "${OUT_DIR}" \
  --min-score 4 \
  --max-ribbons-per-homolog 4

for ext in svg pdf png; do
  [[ -s "${OUT_DIR}/multispecies_A1B1_pcam.pairwise.reference_external.basepair.${ext}" ]] || {
    echo "[ERROR] Missing output: ${OUT_DIR}/multispecies_A1B1_pcam.pairwise.reference_external.basepair.${ext}" >&2
    exit 1
  }
done
echo "[OK] Reference-external base-pair-scaled figure written to ${OUT_DIR}"
