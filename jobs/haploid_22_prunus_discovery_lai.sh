#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# TARGET=all runs targets sequentially; TARGET=<name> supports parallel jobs;
# TARGET=summary aggregates completed per-target manifests.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/prunus_pan_te/discovery}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/analysis/prunus_pan_te/discovery_lai}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
PREP_SCRIPT="${PREP_SCRIPT:-${ROOT_DIR}/scripts/prepare_lai_inputs.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
LAI_ENV="${LAI_ENV:-LAI}"
TARGET="${TARGET:-all}"
THREADS="${LAI_THREADS:-24}"
FORCE_LAI="${FORCE_LAI:-1}"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"
GFF2BED_RESOLVED="${GFF2BED_SCRIPT}"

require_file() {
    local path="$1"
    local label="$2"
    [ -s "${path}" ] || { echo "ERROR: ${label} does not exist or is empty: ${path}" >&2; exit 1; }
}

aggregate_manifests() {
    local summary="${OUT_ROOT}/discovery_lai_summary.tsv"
    local manifest
    shopt -s nullglob
    local manifests=("${OUT_ROOT}"/*/*.lai_manifest.tsv)
    shopt -u nullglob
    [ "${#manifests[@]}" -gt 0 ] || { echo "ERROR: no per-target LAI manifests found under ${OUT_ROOT}" >&2; exit 1; }
    mkdir -p "${OUT_ROOT}"
    printf 'target\tedta_dir\tlai_dir\tgenome\tintact_ltr\tall_ltr\tlai_result\tinput_method\tgenome_bp\twhole_genome_lai\tstatus\n' > "${summary}"
    for manifest in "${manifests[@]}"; do
        tail -n +2 "${manifest}" >> "${summary}"
    done
    local sorted_summary="${summary}.tmp.$$"
    {
        head -n 1 "${summary}"
        tail -n +2 "${summary}" | sort -t $'\t' -k1,1
    } > "${sorted_summary}"
    mv "${sorted_summary}" "${summary}"
    echo "=== LAI summary written: ${summary} ==="
}

run_target() {
    local target="$1"
    local edta_dir="${DISCOVERY_ROOT}/${target}"
    local lai_dir="${OUT_ROOT}/${target}"
    local gff
    local genome intact all result method score genome_bp

    require_file "${PREP_SCRIPT}" "LAI preparation script"
    [ -d "${edta_dir}" ] || { echo "ERROR: discovery directory is missing: ${edta_dir}" >&2; exit 1; }
    gff="$(find "${edta_dir}" -maxdepth 1 -type f \( -name '*.mod.EDTA.TEanno.gff3' -o -name '*.EDTA.TEanno.gff3' \) -size +0c -print -quit)"
    [ -n "${gff}" ] || { echo "ERROR: no EDTA TEanno GFF3 found in ${edta_dir}" >&2; exit 1; }

    mkdir -p "${lai_dir}"
    echo "=== Prepare LAI inputs: ${target} ==="
    python3 "${PREP_SCRIPT}" \
        --edta-dir "${edta_dir}" \
        --out-dir "${lai_dir}" \
        --prefix "${target}" \
        --gff2bed "${GFF2BED_RESOLVED}" \
        --method auto

    require_file "${lai_dir}/${target}.lai_inputs.env" "${target} LAI environment file"
    # shellcheck disable=SC1090
    source "${lai_dir}/${target}.lai_inputs.env"
    genome="${GENOME}"
    intact="${INTACT}"
    all="${ALL}"
    result="${LAI_RESULT}"
    method="${LAI_INPUT_METHOD}"
    require_file "${genome}" "${target} LAI genome"
    require_file "${intact}" "${target} intact LTR pass.list"
    require_file "${all}" "${target} all-LTR input"

    if [ "${FORCE_LAI}" = "1" ] && [ -s "${result}" ]; then
        mv "${result}" "${result}.previous_$(date +%Y%m%d_%H%M%S)_$$"
    fi
    if [ ! -s "${result}" ]; then
        echo "=== Calculate LAI: ${target} ==="
        (
            cd "${lai_dir}"
            LAI -genome "${genome}" -intact "${intact}" -all "${all}" -q -t "${THREADS}"
        )
    else
        echo "Existing LAI output found; skipping: ${result}"
    fi
    require_file "${result}" "${target} LAI result"

    genome_bp="$(awk 'BEGIN{FS="\t"} NR==2 {print $3; exit}' "${result}")"
    score="$(awk 'BEGIN{FS="\t"} $1=="whole_genome" {print $7; exit}' "${result}")"
    genome_bp="${genome_bp:-NA}"
    score="${score:-NA}"
    printf 'target\tedta_dir\tlai_dir\tgenome\tintact_ltr\tall_ltr\tlai_result\tinput_method\tgenome_bp\twhole_genome_lai\tstatus\n' > "${lai_dir}/${target}.lai_manifest.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' \
        "${target}" "${edta_dir}" "${lai_dir}" "${genome}" "${intact}" "${all}" \
        "${result}" "${method}" "${genome_bp}" "${score}" >> "${lai_dir}/${target}.lai_manifest.tsv"
    echo "=== ${target} LAI completed: ${score} ==="
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${TARGET_CONFIG}" "pan-TE target configuration"
source "${TARGET_CONFIG}"

if [ "${TARGET}" = "summary" ]; then
    aggregate_manifests
    exit 0
fi

source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
if [ -z "${GFF2BED_SCRIPT}" ]; then
    for candidate in \
        "${CONDA_PREFIX}/share/EDTA/bin/gff2bed.pl" \
        "${CONDA_PREFIX}/share/LTR_retriever/bin/gff2bed.pl"; do
        if [ -s "${candidate}" ]; then
            GFF2BED_RESOLVED="${candidate}"
            break
        fi
    done
fi
if [ -z "${GFF2BED_RESOLVED}" ] && command -v gff2bed.pl >/dev/null 2>&1; then
    GFF2BED_RESOLVED="$(command -v gff2bed.pl)"
fi
require_file "${GFF2BED_RESOLVED}" "gff2bed.pl"

conda activate "${LAI_ENV}"
command -v LAI >/dev/null 2>&1 || { echo "ERROR: LAI was not found in ${LAI_ENV}" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 was not found in ${LAI_ENV}" >&2; exit 1; }
mkdir -p "${OUT_ROOT}"

if [ "${TARGET}" = "all" ]; then
    for target in "${PAN_TE_TARGETS[@]}"; do
        run_target "${target}"
    done
    aggregate_manifests
else
    target_found=0
    for target in "${PAN_TE_TARGETS[@]}"; do
        if [ "${TARGET}" = "${target}" ]; then
            target_found=1
            run_target "${target}"
            break
        fi
    done
    [ "${target_found}" = "1" ] || { echo "ERROR: unsupported TARGET=${TARGET}" >&2; exit 1; }
fi

echo "=== Discovery LAI job completed ==="
