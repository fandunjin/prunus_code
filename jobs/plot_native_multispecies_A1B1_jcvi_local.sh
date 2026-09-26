#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Render the A1/B1 multi-species ribbon figure with JCVI only.
# Run in WSL from the repository root:
# bash jobs/plot_native_multispecies_A1B1_jcvi_local.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIG_DIR="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_multispecies_A1B1_pcam"
RUNTIME_SHIM="${ROOT}/tmp/jcvi_native_runtime"
JCVI_PYTHON="${JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python3}"

require_file() {
    [[ -s "$1" ]] || { echo "[ERROR] Missing or empty file: $1" >&2; exit 1; }
}

require_file "${FIG_DIR}/multispecies_A1B1.seqids"
require_file "${FIG_DIR}/multispecies_A1B1.layout"
require_file "${JCVI_PYTHON}"

cd "${FIG_DIR}"

# Default output deliberately excludes chromosome circles/labels. This retains
# only native JCVI chromosomes, syntenic ribbons, and necessary row labels.
for format in png pdf; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        multispecies_A1B1.seqids multispecies_A1B1.layout \
        --figsize=18x24 --format="${format}" --notex --nocircles \
        --outfile="multispecies_A1B1_pcam.native.${format}"
done

echo "[OK] Native JCVI outputs:"
echo "  ${FIG_DIR}/multispecies_A1B1_pcam.native.png"
echo "  ${FIG_DIR}/multispecies_A1B1_pcam.native.pdf"
