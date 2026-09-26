#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
AUDIT_SCRIPT="${AUDIT_SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_component_audit.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
MIN_SCORE="${MIN_SCORE:-50}"
MAX_PERIOD="${MAX_PERIOD:-2000}"
FORCE_TRF="${FORCE_TRF:-0}"

[ -n "${TARGET}" ] || { echo "ERROR: set TARGET" >&2; exit 1; }
[ -s "${TARGET_CONFIG}" ] && [ -s "${AUDIT_SCRIPT}" ] && [ -s "${CONDA_SH}" ] || { echo "ERROR: missing script/config" >&2; exit 1; }
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v trf >/dev/null 2>&1 || { echo "ERROR: trf not found in ${EDTA_ENV}" >&2; exit 1; }
set_pan_te_target "${TARGET}"

genome="${OUT_ROOT}/mask/${TARGET}/${TARGET}.fa"
outdir="${FOLLOWUP_ROOT}/tandem/${TARGET}"
raw="${outdir}/${TARGET}.trf.raw"
table="${outdir}/${TARGET}.trf.tsv"
[ -s "${genome}" ] || { echo "ERROR: missing ${genome}" >&2; exit 1; }
mkdir -p "${outdir}"
if [ "${FORCE_TRF}" = "1" ] || [ ! -s "${raw}" ]; then
    trf "${genome}" 2 7 7 80 10 "${MIN_SCORE}" "${MAX_PERIOD}" -ngs -h > "${raw}"
fi
python "${AUDIT_SCRIPT}" trf-to-bed --input "${raw}" --output "${table}" --min-score "${MIN_SCORE}"
[ -s "${table}" ] || { echo "ERROR: empty TRF table ${table}" >&2; exit 1; }
printf "target\tgenome\ttrf_raw\ttrf_table\tmin_score\tmax_period\n%s\t%s\t%s\t%s\t%s\t%s\n" "${TARGET}" "${genome}" "${raw}" "${table}" "${MIN_SCORE}" "${MAX_PERIOD}" > "${outdir}/${TARGET}.tandem_manifest.tsv"
echo "=== Tandem TRF table ready: ${TARGET} ==="
