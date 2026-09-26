#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_ROOT="${OUT_ROOT:-${PAN_ROOT}/figure1_tandem}"
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SUMMARY_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_summary.py"
PARTITION_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_main_summary.py"
MOTIF_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_motif_matrix.py"
MANIFEST="${PAN_ROOT}/summary/standardized_mask_targets.tsv"
TANDEM_ROOT="${PAN_ROOT}/followup/tandem"
SCORE="${SCORE:-50}"

for path in "${PYTHON}" "${SUMMARY_SCRIPT}" "${PARTITION_SCRIPT}" "${MOTIF_SCRIPT}" "${MANIFEST}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done

mkdir -p "${OUT_ROOT}/score${SCORE}" "${OUT_ROOT}/main" "${OUT_ROOT}/main/motifs"
"${PYTHON}" "${SUMMARY_SCRIPT}" \
    --mask-manifest "${MANIFEST}" \
    --tandem-root "${TANDEM_ROOT}" \
    --output-dir "${OUT_ROOT}/score${SCORE}" \
    --window-bp 1000000 \
    --top-motifs-per-target 100 \
    --top-arrays-per-target 100 \
    --min-score "${SCORE}"

"${PYTHON}" "${PARTITION_SCRIPT}" \
    --mask-manifest "${MANIFEST}" \
    --tandem-root "${TANDEM_ROOT}" \
    --tandem-genome-summary "${OUT_ROOT}/score${SCORE}/tandem_genome_summary.tsv" \
    --output-dir "${OUT_ROOT}/main" \
    --min-trf-score "${SCORE}" \
    --reference-group PZI

# This is a complete shared catalogue for short motifs, not the invalid local
# top-N motif list used by the earlier draft heatmap.
"${PYTHON}" "${MOTIF_SCRIPT}" \
    --tandem-root "${TANDEM_ROOT}" \
    --output-dir "${OUT_ROOT}/main/motifs" \
    --min-score "${SCORE}" \
    --max-period-bp 100 \
    --min-total-union-bp 100000

awk -F'\t' 'NR==1 || $9 != 1 || $11 != 1 || $12 != 0 || $13 != 0 {print}' "${OUT_ROOT}/main/tandem_main_partition_qc.tsv"
awk -F'\t' 'NR>1 && ($9 != 1 || $11 != 1 || $12 != 0 || $13 != 0) {bad=1} END {exit bad}' "${OUT_ROOT}/main/tandem_main_partition_qc.tsv"
echo "=== Main-figure tandem partition completed ==="
