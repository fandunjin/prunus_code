#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
TARGET="${TARGET:-}"
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_loci.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-minimap}"

case "${TARGET}" in
    pziA1) HAP=hapA1 ;;
    pziA2) HAP=hapA2 ;;
    pziB1) HAP=hapB1 ;;
    pziB2) HAP=hapB2 ;;
    *) echo "ERROR: TARGET must be pziA1, pziA2, pziB1 or pziB2" >&2; exit 1 ;;
esac
BAM="${BAM:-${ROOT_DIR}/analysis/pzi_expansion_qc/${HAP}/ont_depth/pzi_${HAP}.ont.primary.mapq20.sorted.bam}"
for path in "${PYTHON}" "${SCRIPT}" "${FIG_ROOT}/loci/pzi_tandem_array_candidates.tsv" "${BAM}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${SAMTOOLS_ENV}"
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools was not available" >&2; exit 1; }
mkdir -p "${FIG_ROOT}/ont"
"${PYTHON}" "${SCRIPT}" ont \
    --candidates "${FIG_ROOT}/loci/pzi_tandem_array_candidates.tsv" \
    --target "${TARGET}" \
    --bam "${BAM}" \
    --output "${FIG_ROOT}/ont/${TARGET}.ont_support.tsv" \
    --top-arrays 15 \
    --flank-bp 1000 \
    --min-mapq 20 \
    --min-junction-reads 3
echo "=== Tandem-array ONT junction support completed: ${TARGET} ==="
