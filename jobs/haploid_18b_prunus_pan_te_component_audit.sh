#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
POSTPROCESS_DIR="${POSTPROCESS_DIR:-${OUT_ROOT}/postprocess_v2}"
PAN_LIBRARY="${PAN_LIBRARY:-${OUT_ROOT}/library/frozen/prunus_pan_te_v1.fa}"
OUTDIR="${OUTDIR:-${FOLLOWUP_ROOT}/component_audit}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_component_audit.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
THREADS="${THREADS:-4}"
MAX_COMPONENTS="${MAX_COMPONENTS:-50}"
CUMULATIVE_THRESHOLD="${CUMULATIVE_THRESHOLD:-0.60}"

for path in "${SCRIPT}" "${PAN_LIBRARY}" "${POSTPROCESS_DIR}/te_driver_families.tsv" "${CONDA_SH}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v TEsorter >/dev/null 2>&1 || { echo "ERROR: TEsorter was not found" >&2; exit 1; }
mkdir -p "${OUTDIR}"

python "${SCRIPT}" select \
    --driver-table "${POSTPROCESS_DIR}/te_driver_families.tsv" \
    --pan-library "${PAN_LIBRARY}" \
    --output-dir "${OUTDIR}" \
    --max-components "${MAX_COMPONENTS}" \
    --cumulative-threshold "${CUMULATIVE_THRESHOLD}"

selected_fasta="${OUTDIR}/selected_components.fa"
tesorter_cls="${selected_fasta}.rexdb-plant.cls.tsv"
if [ ! -s "${tesorter_cls}" ]; then
    TEsorter "${selected_fasta}" -db rexdb-plant -p "${THREADS}"
fi

python "${SCRIPT}" audit \
    --selected-components "${OUTDIR}/selected_components.tsv" \
    --mask-root "${OUT_ROOT}/mask" \
    --tandem-root "${FOLLOWUP_ROOT}/tandem" \
    --tesorter-cls "${tesorter_cls}" \
    --output-dir "${OUTDIR}"

echo "=== Component audit completed ==="
echo "Review and edit: ${OUTDIR}/component_decisions.tsv"
