#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
TARGET="${TARGET:-}"
THREADS="${THREADS:-8}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
FLANKS="${FIG_ROOT}/loci/pzi_tandem_array_flanks.fa"
OUT="${FIG_ROOT}/loci/paf"
FORCE_RERUN="${FORCE_RERUN:-0}"

case "${TARGET}" in
    pziA1|pziA2|pziB1|pziB2|pspA1|pspA2|pspB1|pspB2|ppaA|ppaB|pcam|ppe) ;;
    *) echo "ERROR: TARGET is not a standardized genomic package" >&2; exit 1 ;;
esac
for path in "${CONDA_SH}" "${FLANKS}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
TARGET_FASTA="${PAN_ROOT}/mask/${TARGET}/${TARGET}.fa"
[ -s "${TARGET_FASTA}" ] || { echo "ERROR: missing standardized target FASTA ${TARGET_FASTA}" >&2; exit 1; }
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 was not available" >&2; exit 1; }
mkdir -p "${OUT}"
PAF="${OUT}/all_pzi_array_flanks_to_${TARGET}.paf"
[ "${FORCE_RERUN}" = 1 ] && rm -f "${PAF}"
if [ ! -s "${PAF}" ]; then
    minimap2 -x asm20 -N 5 -t "${THREADS}" "${TARGET_FASTA}" "${FLANKS}" > "${PAF}"
fi
[ -s "${PAF}" ] || { echo "ERROR: empty PAF ${PAF}" >&2; exit 1; }
printf 'target\ttarget_fasta\tflanks\tpaf\n%s\t%s\t%s\t%s\n' "${TARGET}" "${TARGET_FASTA}" "${FLANKS}" "${PAF}" > "${OUT}/all_pzi_array_flanks_to_${TARGET}.manifest.tsv"
echo "=== Tandem flank mapping completed: ${TARGET} ==="
