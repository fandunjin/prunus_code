#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare and render the 0903 seven-species ribbon figure with native JCVI.
# Run from WSL at the repository root.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREPARE="${ROOT}/scripts/plot_0903_multispecies_A1B1_ribbon.py"
JCVI_PYTHON="${JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python3}"
RUNTIME_SHIM="${ROOT}/tmp/jcvi_native_runtime"
# Reuse the reference directory so the 0903 native outputs replace the
# previous approximation under the same user-facing filenames.
OUT_DIR="${ROOT}/analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_multispecies_A1B1_pcam_pairwise"

require_file() { [[ -s "$1" ]] || { echo "[ERROR] Missing or empty file: $1" >&2; exit 1; }; }
require_file "${PREPARE}"
require_file "${JCVI_PYTHON}"

# Generate the 0903-derived JCVI input bundle.  The renderer below is the
# native JCVI command; the preparation script only stages BED/anchor metadata.
"${JCVI_PYTHON}" "${PREPARE}" \
    --hap-input "${ROOT}/0903/haploid_repartition/inputs" \
    --jcvi-input "${ROOT}/0903/downstream_analysis/figure2_repartition_ks/inputs" \
    --pair-work "${ROOT}/0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work" \
    --rch-input "${ROOT}/downstream_analysis/local_A1B1_ribbon_inputs" \
    --outdir "${OUT_DIR}" \
    --min-score 4 \
    --max-ribbons-per-homolog 4 \
    --prepare-only

# Guard against accidentally rendering a base-pair figure with the ordinary
# gene-normalized layout.  The two layouts must differ in their xstart/xend
# track spans, and the scale table records the source BED spans used.
if cmp -s "${OUT_DIR}/multispecies_0903_A1B1_pcam.layout" \
        "${OUT_DIR}/multispecies_0903_A1B1_pcam_basepair.layout"; then
    echo "[ERROR] Base-pair layout is identical to the normalized layout" >&2
    exit 1
fi
require_file "${OUT_DIR}/basepair_track_scale.tsv"

cd "${OUT_DIR}"

for format in png pdf svg; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        multispecies_0903_A1B1_pcam.seqids multispecies_0903_A1B1_pcam.layout \
        --figsize=18x24 --format="${format}" --notex --keep-chrlabels \
        --chrstyle=roundrect --shadestyle=curve \
        --outfile="multispecies_A1B1_pcam.pairwise.native.${format}"
done

for format in png pdf svg; do
    PYTHONPATH="${RUNTIME_SHIM}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${JCVI_PYTHON}" -m jcvi.graphics.karyotype \
        multispecies_0903_A1B1_pcam.seqids multispecies_0903_A1B1_pcam_basepair.layout \
        --figsize=18x24 --format="${format}" --notex --nocircles --basepair --keep-chrlabels \
        --chrstyle=roundrect --shadestyle=curve \
        --outfile="multispecies_A1B1_pcam.pairwise.native.basepair.${format}"
done

for format in png pdf svg; do
    require_file "${OUT_DIR}/multispecies_A1B1_pcam.pairwise.native.${format}"
    require_file "${OUT_DIR}/multispecies_A1B1_pcam.pairwise.native.basepair.${format}"
done

echo "[OK] Native JCVI figure: ${OUT_DIR}/multispecies_A1B1_pcam.pairwise.native.pdf"
echo "[OK] Native JCVI base-pair figure: ${OUT_DIR}/multispecies_A1B1_pcam.pairwise.native.basepair.pdf"
