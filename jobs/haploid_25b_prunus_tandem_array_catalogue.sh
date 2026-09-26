#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_loci.py"
MANIFEST="${PAN_ROOT}/summary/standardized_mask_targets.tsv"
TANDEM_ROOT="${PAN_ROOT}/followup/tandem"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-haphic}"

for path in "${PYTHON}" "${SCRIPT}" "${MANIFEST}" "${CONDA_SH}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${CONDA_SH}"
conda activate "${SAMTOOLS_ENV}"
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools was not available" >&2; exit 1; }
for target in pziA1 pziA2 pziB1 pziB2; do
    fasta="${PAN_ROOT}/mask/${target}/${target}.fa"
    [ -s "${fasta}" ] || { echo "ERROR: missing ${fasta}" >&2; exit 1; }
    [ -s "${fasta}.fai" ] || samtools faidx "${fasta}"
done

"${PYTHON}" "${SCRIPT}" catalogue \
    --mask-manifest "${MANIFEST}" \
    --tandem-root "${TANDEM_ROOT}" \
    --output-dir "${FIG_ROOT}/loci" \
    --min-score 50 \
    --join-gap-bp 10000 \
    --min-array-union-bp 100000 \
    --top-arrays-per-target 25 \
    --flank-bp 5000

echo "=== PZI tandem-array catalogue completed ==="
