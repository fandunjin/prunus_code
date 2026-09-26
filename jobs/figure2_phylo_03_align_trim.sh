#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
SOURCE="${BASE}/03_single_copy"
INPUTS="${BASE}/01_frozen_inputs"
OUT="${BASE}/04_alignments"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-48}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

test -e "${SOURCE}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_02 first." >&2; exit 3; }
test ! -e "${OUT}" || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 3; }
mkdir -p "${OUT}"
set +u
source "${CONDA_SH}"
conda activate "${TREE_ENV}"
set -u
for exe in mafft trimal parallel; do command -v "${exe}" >/dev/null || { echo "[ERROR] ${exe} unavailable" >&2; exit 2; }; done

export MAFFT_BIN="$(command -v mafft)"
export TRIMAL_BIN="$(command -v trimal)"
for dataset in main sensitivity; do
    mkdir -p "${OUT}/${dataset}/mafft" "${OUT}/${dataset}/trimmed" "${OUT}/${dataset}/logs"
    find "${SOURCE}/${dataset}/unaligned" -maxdepth 1 -type f -name '*.fa' -print0 | sort -z | \
        parallel -0 -j "${THREADS}" --halt now,fail=1 \
        'base=$(basename {} .fa); "${MAFFT_BIN}" --localpair --maxiterate 1000 --thread 1 {} > "'"${OUT}/${dataset}/mafft"'/${base}.fa" 2> "'"${OUT}/${dataset}/logs"'/${base}.mafft.log"; "${TRIMAL_BIN}" -in "'"${OUT}/${dataset}/mafft"'/${base}.fa" -out "'"${OUT}/${dataset}/trimmed"'/${base}.fa" -fasta -automated1 > "'"${OUT}/${dataset}/logs"'/${base}.trimal.log" 2>&1'

    find "${INPUTS}/${dataset}/pep" -maxdepth 1 -type f -name '*.fa' -printf '%f\n' | sed 's/\.fa$//' | sort > "${OUT}/${dataset}/species.txt"
    "${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_alignment_qc.py" \
        --alignment-dir "${OUT}/${dataset}/trimmed" \
        --species-file "${OUT}/${dataset}/species.txt" \
        --out-dir "${OUT}/${dataset}/qc" \
        --min-length 100 \
        --min-informative 10 \
        --max-missing 0.20 \
        --min-alignments 200
done
{
    mafft --version
    trimal --version
    parallel --version | head -n 1
} > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
