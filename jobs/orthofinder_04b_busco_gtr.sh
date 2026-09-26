#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/orthofinder/species/pep}"
OUT_DIR="${OUT_DIR:-${PEP_DIR}/busco_eudicots}"
SAFE_DIR="${SAFE_DIR:-${PEP_DIR}/busco_safe_inputs}"
SUMMARY="${SUMMARY:-${PEP_DIR}/busco_eudicots_summary.tsv}"
INPUT_FASTA="${INPUT_FASTA:-${ROOT_DIR}/data/gtr/Gtr.NCBI.prot.faa}"
LABEL="${LABEL:-gtr}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
FORCE="${FORCE:-0}"
SAFE_FASTA_SCRIPT="${SAFE_FASTA_SCRIPT:-${ROOT_DIR}/scripts/downstream_make_busco_safe_fasta.py}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${INPUT_FASTA}" "GTR protein FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
require_file "${SAFE_FASTA_SCRIPT}" "BUSCO-safe FASTA helper"
mkdir -p "${OUT_DIR}" "${SAFE_DIR}" "${ROOT_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${BUSCO_ENV}"
set -u

command -v busco >/dev/null 2>&1 || {
    echo "ERROR: busco was not found in BUSCO_ENV=${BUSCO_ENV}." >&2
    exit 1
}

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN=python3
else
    PYTHON_BIN=python
fi

run_name="${LABEL}_protein_eudicots"
run_dir="${OUT_DIR}/${run_name}"
safe_pep="${SAFE_DIR}/${LABEL}.busco_safe.pep.fa"
safe_map="${SAFE_DIR}/${LABEL}.busco_safe_header_map.tsv"

echo "=== BUSCO protein-mode evaluation: ${LABEL} ==="
echo "Input protein FASTA : ${INPUT_FASTA}"
echo "BUSCO output        : ${OUT_DIR}"
echo "Lineage             : ${BUSCO_LINEAGE}"
echo "Threads             : ${THREADS}"
echo "Force               : ${FORCE}"
echo

summary_file=$(find "${run_dir}" -maxdepth 3 -name 'short_summary*.txt' -type f -size +0c 2>/dev/null | head -n 1 || true)
if [ -n "${summary_file}" ] && [ "${FORCE}" != "1" ]; then
    echo "Skip ${LABEL}: existing summary ${summary_file}"
else
    "${PYTHON_BIN}" "${SAFE_FASTA_SCRIPT}" \
        --input "${INPUT_FASTA}" \
        --output "${safe_pep}" \
        --map "${safe_map}" \
        --label "${LABEL}"

    echo ">>> BUSCO ${LABEL}"
    busco_args=(
        -i "${safe_pep}"
        -m proteins
        -l "${BUSCO_LINEAGE}"
        -c "${THREADS}"
        --offline
        -o "${run_name}"
        --out_path "${OUT_DIR}"
    )
    if [ "${FORCE}" = "1" ]; then
        busco_args+=(--force)
    elif [ -d "${run_dir}" ]; then
        echo "Resume incomplete BUSCO run: ${run_dir}"
        busco_args+=(--restart)
    fi
    busco "${busco_args[@]}"
fi

summary_file=$(find "${run_dir}" -maxdepth 3 -name 'short_summary*.txt' -type f -size +0c 2>/dev/null | head -n 1 || true)
if [ -z "${summary_file}" ]; then
    echo "ERROR: no BUSCO short_summary found under ${run_dir}" >&2
    exit 1
fi
busco_line=$(grep -m1 -E 'C:[0-9.]+%' "${summary_file}" || true)
if [ -z "${busco_line}" ]; then
    echo "ERROR: could not parse BUSCO summary: ${summary_file}" >&2
    exit 1
fi
complete=$(printf '%s\n' "${busco_line}" | sed -E 's/.*C:([0-9.]+)%.*/\1/')
single=$(printf '%s\n' "${busco_line}" | sed -E 's/.*S:([0-9.]+)%.*/\1/')
duplicated=$(printf '%s\n' "${busco_line}" | sed -E 's/.*D:([0-9.]+)%.*/\1/')
fragmented=$(printf '%s\n' "${busco_line}" | sed -E 's/.*F:([0-9.]+)%.*/\1/')
missing=$(printf '%s\n' "${busco_line}" | sed -E 's/.*M:([0-9.]+)%.*/\1/')

row=$(printf '%s\tPASS\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${LABEL}" "${complete}" "${single}" "${duplicated}" "${fragmented}" "${missing}" "${summary_file}")

if [ ! -s "${SUMMARY}" ]; then
    printf 'label\tstatus\tcomplete_pct\tsingle_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\tbusco_summary\n' > "${SUMMARY}"
fi
if grep -q "^${LABEL}"$'\t' "${SUMMARY}"; then
    tmp_summary="$(mktemp)"
    awk -v r="${row}" -v l="${LABEL}" 'BEGIN{FS=OFS="\t"} $1==l{print r; next} {print}' "${SUMMARY}" > "${tmp_summary}"
    mv "${tmp_summary}" "${SUMMARY}"
    echo "Updated existing ${LABEL} row in ${SUMMARY}"
else
    printf '%s\n' "${row}" >> "${SUMMARY}"
    echo "Appended ${LABEL} row to ${SUMMARY}"
fi

echo
echo "=== BUSCO evaluation completed ==="
echo "Summary: ${SUMMARY}"
echo "New row:"
echo "${row}"
