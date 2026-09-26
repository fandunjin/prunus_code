#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_tandem_repeatmasker_overlap.py}"
OUT="${OUT:-${PAN_ROOT}/followup/tandem_summary/tandem_repeatmasker_overlap.tsv}"
for path in "${PYTHON}" "${SCRIPT}"; do [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }; done
"${PYTHON}" "${SCRIPT}" --tandem-root "${PAN_ROOT}/followup/tandem" --mask-root "${PAN_ROOT}/mask" --output "${OUT}"
echo "=== Tandem RepeatMasker-overlap diagnosis completed ==="
