#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON_ENV="${PYTHON_ENV:-haphic}"
PYTHON_SCRIPT="${ROOT_DIR}/scripts/merge_structure_qc_summaries.py"
OUT="${OUT:-${ROOT_DIR}/annotation/summary/gene_structure_metrics_all.tsv}"
PZI_SUMMARY="${PZI_SUMMARY:-${ROOT_DIR}/annotation/genes/pzi_eviann/postprocess/structure_qc/gene_structure_metrics.tsv}"
PSP_SUMMARY="${PSP_SUMMARY:-${ROOT_DIR}/data/psp/annotation/genes/psp32_eviann/postprocess/structure_qc/gene_structure_metrics.tsv}"
PPA_SUMMARY="${PPA_SUMMARY:-${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/structure_qc/gene_structure_metrics.tsv}"
[[ -s "$CONDA_SH" && -s "$PYTHON_SCRIPT" ]] || { echo "ERROR: conda or merge script missing" >&2; exit 1; }
source "$CONDA_SH"
conda activate "$PYTHON_ENV"
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 unavailable" >&2; exit 1; }
python3 "$PYTHON_SCRIPT" --input "pzi=$PZI_SUMMARY" --input "psp=$PSP_SUMMARY" --input "ppa=$PPA_SUMMARY" --output "$OUT"
