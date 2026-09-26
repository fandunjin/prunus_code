#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_tandem_summary.py}"
MANIFEST="${MANIFEST:-${PAN_ROOT}/summary/standardized_mask_targets.tsv}"
TANDEM_ROOT="${TANDEM_ROOT:-${PAN_ROOT}/followup/tandem}"
OUT="${OUT:-${PAN_ROOT}/followup/tandem_summary}"
for path in "${PYTHON}" "${SCRIPT}" "${MANIFEST}"; do [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }; done
"${PYTHON}" "${SCRIPT}" --mask-manifest "${MANIFEST}" --tandem-root "${TANDEM_ROOT}" --output-dir "${OUT}" --window-bp 1000000 --top-motifs-per-target 100 --top-arrays-per-target 100
echo "=== Tandem-repeat summary completed ==="
