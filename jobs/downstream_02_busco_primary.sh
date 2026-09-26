#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OUT_DIR="${ANALYSIS_DIR}/qc/busco_primary"
SAFE_DIR="${OUT_DIR}/safe_inputs"

mkdir -p "${OUT_DIR}" "${SAFE_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${BUSCO_ENV}"
set -u

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN=python
    else
        echo "[ERROR] Neither ${PYTHON_BIN} nor python is available in BUSCO_ENV=${BUSCO_ENV}" >&2
        exit 2
    fi
fi

tail -n +2 "${ANALYSIS_DIR}/canonical_manifest.tsv" | while IFS=$'\t' read -r label role tier data_class isoform_strategy include_main include_sensitivity include_synteny canonical_dir genome_records genome_bp pep_records pep_aa cds_records cds_bp selected skipped; do
    pep="${canonical_dir}/primary.pep.fa"
    if [ ! -s "${pep}" ]; then
        echo "Skip BUSCO for ${label}: missing primary.pep.fa"
        continue
    fi

    run_name="${label}_primary_protein_eudicots"
    run_dir="${OUT_DIR}/${run_name}"
    if find "${run_dir}" -maxdepth 2 -name 'short_summary*.txt' -size +0c 2>/dev/null | grep -q .; then
        echo "Skip BUSCO for ${label}: existing complete BUSCO summary in ${run_dir}"
        continue
    fi

    safe_pep="${SAFE_DIR}/${label}.busco_safe.pep.fa"
    safe_map="${SAFE_DIR}/${label}.busco_safe_header_map.tsv"
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_make_busco_safe_fasta.py" \
        --input "${pep}" \
        --output "${safe_pep}" \
        --map "${safe_map}" \
        --label "${label}"

    busco \
        -i "${safe_pep}" \
        -m proteins \
        -l "${BUSCO_LINEAGE}" \
        -c "${THREADS}" \
        -o "${run_name}" \
        --out_path "${OUT_DIR}" \
        --force \
        --offline
done
