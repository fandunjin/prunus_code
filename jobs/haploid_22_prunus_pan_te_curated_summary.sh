#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
DECISIONS="${DECISIONS:-${PAN_ROOT}/followup/component_audit/component_decisions.tsv}"
[ -s "${DECISIONS}" ] || DECISIONS="${PAN_ROOT}/component_decisions.tsv"
DRIVERS="${DRIVERS:-${PAN_ROOT}/postprocess_v2/te_driver_families.tsv}"
INTACT="${INTACT:-${PAN_ROOT}/postprocess_v2/te_intact_ltr.tsv}"
OUT="${OUT:-${PAN_ROOT}/followup/curated_summary}"
[ -s "${DRIVERS}" ] && [ -s "${DECISIONS}" ] && [ -s "${INTACT}" ] || { echo "ERROR: missing curated-summary input" >&2; exit 1; }
"${PYTHON}" "${ROOT_DIR}/scripts/prunus_pan_te_curated_summary.py" --drivers "${DRIVERS}" --decisions "${DECISIONS}" --intact "${INTACT}" --output-dir "${OUT}"
echo "=== curated TE summary completed ==="
