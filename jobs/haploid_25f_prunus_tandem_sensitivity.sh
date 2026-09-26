#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SUMMARY_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_summary.py"
PARTITION_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_main_summary.py"
SENSITIVITY_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_sensitivity.py"
MANIFEST="${PAN_ROOT}/summary/standardized_mask_targets.tsv"
TANDEM_ROOT="${PAN_ROOT}/followup/tandem"

for path in "${PYTHON}" "${SUMMARY_SCRIPT}" "${PARTITION_SCRIPT}" "${SENSITIVITY_SCRIPT}" "${MANIFEST}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
mkdir -p "${FIG_ROOT}/sensitivity"
args=()
for score in 50 80 120; do
    out="${FIG_ROOT}/sensitivity/score${score}"
    mkdir -p "${out}"
    "${PYTHON}" "${SUMMARY_SCRIPT}" --mask-manifest "${MANIFEST}" --tandem-root "${TANDEM_ROOT}" --output-dir "${out}" --window-bp 1000000 --top-motifs-per-target 100 --top-arrays-per-target 100 --min-score "${score}"
    "${PYTHON}" "${PARTITION_SCRIPT}" --mask-manifest "${MANIFEST}" --tandem-root "${TANDEM_ROOT}" --tandem-genome-summary "${out}/tandem_genome_summary.tsv" --output-dir "${out}/partition" --min-trf-score "${score}" --reference-group PZI
    args+=(--input "${score}:${out}/partition/tandem_main_partition.tsv")
done
"${PYTHON}" "${SENSITIVITY_SCRIPT}" "${args[@]}" --output "${FIG_ROOT}/main/tandem_score_sensitivity.tsv"
echo "=== Tandem score sensitivity completed ==="
