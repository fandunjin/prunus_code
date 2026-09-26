#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${OUT_DIR:-${PEP_DIR}/busco_eudicots}"
SAFE_DIR="${SAFE_DIR:-${PEP_DIR}/busco_safe_inputs}"
SUMMARY="${SUMMARY:-${PEP_DIR}/busco_eudicots_summary.tsv}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
FORCE="${FORCE:-0}"
SAFE_FASTA_SCRIPT="${ROOT_DIR}/scripts/downstream_make_busco_safe_fasta.py"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
require_file "${SAFE_FASTA_SCRIPT}" "BUSCO-safe FASTA helper"
mkdir -p "${OUT_DIR}" "${SAFE_DIR}" "${ROOT_DIR}/logs"

mapfile -t PEP_FILES < <(find "${PEP_DIR}" -maxdepth 1 -type f -name '*.pep.fa' -size +0c | sort)
if [ "${#PEP_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no non-empty *.pep.fa files found in ${PEP_DIR}" >&2
    exit 1
fi

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

echo "=== BUSCO protein-mode evaluation for OrthoFinder species set ==="
echo "Input directory : ${PEP_DIR}"
echo "Protein files   : ${#PEP_FILES[@]}"
echo "BUSCO output    : ${OUT_DIR}"
echo "Safe inputs     : ${SAFE_DIR}"
echo "Lineage         : ${BUSCO_LINEAGE}"
echo "Threads         : ${THREADS}"
echo "Force           : ${FORCE}"
echo

for pep in "${PEP_FILES[@]}"; do
    label="$(basename "${pep}" .pep.fa)"
    run_name="${label}_protein_eudicots"
    run_dir="${OUT_DIR}/${run_name}"
    summary_file=$(find "${run_dir}" -maxdepth 3 -name 'short_summary*.txt' -type f -size +0c 2>/dev/null | head -n 1 || true)

    if [ -n "${summary_file}" ] && [ "${FORCE}" != "1" ]; then
        echo "Skip ${label}: existing summary ${summary_file}"
        continue
    fi

    safe_pep="${SAFE_DIR}/${label}.busco_safe.pep.fa"
    safe_map="${SAFE_DIR}/${label}.busco_safe_header_map.tsv"
    "${PYTHON_BIN}" "${SAFE_FASTA_SCRIPT}" \
        --input "${pep}" \
        --output "${safe_pep}" \
        --map "${safe_map}" \
        --label "${label}"

    echo ">>> BUSCO ${label}"
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
done

printf 'label\tstatus\tcomplete_pct\tsingle_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\tbusco_summary\n' > "${SUMMARY}"
for pep in "${PEP_FILES[@]}"; do
    label="$(basename "${pep}" .pep.fa)"
    run_dir="${OUT_DIR}/${label}_protein_eudicots"
    summary_file=$(find "${run_dir}" -maxdepth 3 -name 'short_summary*.txt' -type f -size +0c 2>/dev/null | head -n 1 || true)
    if [ -z "${summary_file}" ]; then
        printf '%s\tFAIL\tNA\tNA\tNA\tNA\tNA\tNA\n' "${label}" >> "${SUMMARY}"
        continue
    fi
    busco_line=$(grep -m1 -E 'C:[0-9.]+%' "${summary_file}" || true)
    if [ -z "${busco_line}" ]; then
        printf '%s\tFAIL\tNA\tNA\tNA\tNA\tNA\t%s\n' "${label}" "${summary_file}" >> "${SUMMARY}"
        continue
    fi
    complete=$(printf '%s\n' "${busco_line}" | sed -E 's/.*C:([0-9.]+)%.*/\1/')
    single=$(printf '%s\n' "${busco_line}" | sed -E 's/.*S:([0-9.]+)%.*/\1/')
    duplicated=$(printf '%s\n' "${busco_line}" | sed -E 's/.*D:([0-9.]+)%.*/\1/')
    fragmented=$(printf '%s\n' "${busco_line}" | sed -E 's/.*F:([0-9.]+)%.*/\1/')
    missing=$(printf '%s\n' "${busco_line}" | sed -E 's/.*M:([0-9.]+)%.*/\1/')
    printf '%s\tPASS\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${label}" "${complete}" "${single}" "${duplicated}" "${fragmented}" "${missing}" "${summary_file}" >> "${SUMMARY}"
done

echo
echo "=== BUSCO evaluation completed ==="
echo "Summary: ${SUMMARY}"
cat "${SUMMARY}"
