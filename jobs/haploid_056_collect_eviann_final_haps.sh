#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PY_ENV="${PY_ENV:-haphic}"
EVIANN_ROOT="${ROOT_DIR}/annotation/genes/final_haplotypes_eviann"
EVIANN_PACKAGE="${ROOT_DIR}/final_haplotype_genomes_eviann"
COLLECT_SCRIPT="${ROOT_DIR}/scripts/collect_eviann_haplotype_genomes.py"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${COLLECT_SCRIPT}" "EviAnn final haplotype collection script"

for hap in hapA1 hapA2 hapB1 hapB2; do
    require_file "${ROOT_DIR}/haploid/final_haplotypes/pzi_${hap}.fa" "${hap} genome FASTA"
    require_file "${EVIANN_ROOT}/${hap}/postprocess/pzi_${hap}.eviann.genes.gff3" "${hap} postprocessed GFF3"
    require_file "${EVIANN_ROOT}/${hap}/postprocess/pzi_${hap}.eviann.cds.fa" "${hap} postprocessed CDS"
    require_file "${EVIANN_ROOT}/${hap}/postprocess/pzi_${hap}.eviann.protein.fa" "${hap} postprocessed protein"
done

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Collect EviAnn-refined final haplotypes ==="
echo "EviAnn root : ${EVIANN_ROOT}"
echo "Package dir : ${EVIANN_PACKAGE}"
echo

"${PYTHON_BIN}" "${COLLECT_SCRIPT}" \
    --root-dir "${ROOT_DIR}" \
    --eviann-root "${EVIANN_ROOT}" \
    --out-dir "${EVIANN_PACKAGE}"

echo
echo "=== Collection finished ==="
echo "${EVIANN_PACKAGE}"
