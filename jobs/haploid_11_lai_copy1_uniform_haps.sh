#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PAN_RESULT_ROOT="${ROOT_DIR}/annotation/repeats/final_haplotypes/panEDTA_uniform"
OUT_ROOT="${ROOT_DIR}/annotation/repeats/lai_copy1_uniform_haps"
LOG_DIR="${ROOT_DIR}/logs"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
LAI_ENV="${LAI_ENV:-LAI}"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_lai_inputs.py"
THREADS="${THREADS:-24}"
FORCE_LAI="${FORCE_LAI:-0}"
COPY1_MONO="${COPY1_MONO:-1}"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"

HAPS=(copy1 hapA1 hapA2 hapB1 hapB2)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

stage_uniform_lai_inputs() {
    local label="$1"
    local source_dir="$2"
    local base="$3"
    local prefix="$4"
    local stage_dir="${OUT_ROOT}/uniform_edta_inputs/${label}"
    local gff_source="${source_dir}/${base}.mod.EDTA.TEanno.gff3"
    local genome_source="${source_dir}/${base}.mod"
    local rmout_source="${source_dir}/${base}.mod.panEDTA.out"
    local raw_dir="${source_dir}/${base}.mod.EDTA.raw"
    local pass_source=""

    require_file "${gff_source}" "${label} unified EDTA GFF3"
    require_file "${genome_source}" "${label} EDTA-normalized genome"
    require_file "${rmout_source}" "${label} unified panEDTA RepeatMasker output"
    mkdir -p "${stage_dir}"
    ln -sfn "${gff_source}" "${stage_dir}/${base}.mod.EDTA.TEanno.gff3"
    ln -sfn "${genome_source}" "${stage_dir}/${base}.mod"

    if [[ -d "${raw_dir}" ]]; then
        pass_source="$(find "${raw_dir}" -type f -name '*.real.pass.list' -size +0c -print -quit 2>/dev/null || true)"
    fi
    if [[ -n "${pass_source}" ]]; then
        ln -sfn "${pass_source}" "${stage_dir}/${prefix}.real.pass.list"
        ln -sfn "${rmout_source}" "${stage_dir}/${prefix}.all_ltr.out"
        PREP_METHOD="native"
    else
        PREP_METHOD="edta-gff3"
    fi

    STAGED_EDTA_DIR="${stage_dir}"
    UNIFORM_RMOUT="${rmout_source}"
}

mkdir -p "${OUT_ROOT}" "${LOG_DIR}"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${PREP_SCRIPT}" "LAI input preparation script"
require_file "${COPY1_UNIFORM_EDTA_DIR}/pzi_copy1.fa.mod.EDTA.TEanno.gff3" "copy1 unified EDTA GFF3"

source "${CONDA_SH}"
conda activate "${LAI_ENV}"
command -v LAI >/dev/null 2>&1 || die "LAI is unavailable in ${LAI_ENV}"

if [[ -z "${GFF2BED_SCRIPT}" ]]; then
    if [[ -s "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl" ]]; then
        GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl"
    elif [[ -s "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl" ]]; then
        GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl"
    elif command -v gff2bed.pl >/dev/null 2>&1; then
        GFF2BED_RESOLVED="$(command -v gff2bed.pl)"
    else
        die "gff2bed.pl was not found; set GFF2BED_SCRIPT=/path/to/gff2bed.pl"
    fi
else
    GFF2BED_RESOLVED="${GFF2BED_SCRIPT}"
fi
require_file "${GFF2BED_RESOLVED}" "gff2bed.pl"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

SUMMARY="${OUT_ROOT}/copy1_uniform_haplotypes_lai_summary.tsv"
printf 'label\tedta_dir\tlai_dir\tgenome\tintact_ltr\tall_ltr\tlai_result\tinput_method\tlai_mode\tstatus\n' > "${SUMMARY}"

echo "=== LAI for copy1 and panEDTA-uniform haplotypes ==="
echo "Output      : ${OUT_ROOT}"
echo "LAI env     : ${LAI_ENV}"
echo "Threads     : ${THREADS}"
echo "Force LAI   : ${FORCE_LAI}"
echo "copy1 -mono : ${COPY1_MONO}"
echo "gff2bed     : ${GFF2BED_RESOLVED}"
echo

for label in "${HAPS[@]}"; do
    if [[ "${label}" == "copy1" ]]; then
        prefix="pzi_copy1"
        source_dir="${COPY1_UNIFORM_EDTA_DIR}"
        base="pzi_copy1.fa"
    else
        prefix="pzi_${label}"
        source_dir="${PAN_RESULT_ROOT}"
        base="pzi_${label}.fa"
    fi

    stage_uniform_lai_inputs "${label}" "${source_dir}" "${base}" "${prefix}"
    edta_dir="${STAGED_EDTA_DIR}"
    lai_dir="${OUT_ROOT}/${label}_lai"
    mkdir -p "${lai_dir}"

    echo "=== Prepare LAI inputs: ${label} ==="
    echo "EDTA dir: ${edta_dir}"
    echo "LAI dir : ${lai_dir}"
    "${PYTHON_BIN}" "${PREP_SCRIPT}" \
        --edta-dir "${edta_dir}" \
        --out-dir "${lai_dir}" \
        --prefix "${prefix}" \
        --gff2bed "${GFF2BED_RESOLVED}" \
        --method "${PREP_METHOD}"

    require_file "${lai_dir}/${prefix}.lai_inputs.env" "${label} LAI environment file"
    # shellcheck disable=SC1090
    source "${lai_dir}/${prefix}.lai_inputs.env"
    ALL="${UNIFORM_RMOUT}"
    LAI_RESULT="${lai_dir}/$(basename "${ALL}").LAI"
    LAI_INPUT_METHOD="${LAI_INPUT_METHOD}+panEDTA-rmout"
    require_file "${GENOME}" "${label} LAI genome"
    require_file "${INTACT}" "${label} intact LTR pass.list"
    require_file "${ALL}" "${label} all-LTR RepeatMasker output"

    cd "${lai_dir}"
    LAI_MODE="standard"
    LAI_EXTRA_ARGS=()
    if [[ "${label}" == "copy1" && "${COPY1_MONO}" == "1" ]]; then
        LAI_MODE="mono"
        LAI_EXTRA_ARGS+=(-mono)
    fi
    METHOD_MARKER="${lai_dir}/.${prefix}.uniform_panedta_rmout_lai_${LAI_MODE}"
    if [[ ( "${FORCE_LAI}" == "1" || ! -e "${METHOD_MARKER}" ) && -s "${LAI_RESULT}" ]]; then
        backup="${LAI_RESULT}.previous_$(date +%Y%m%d_%H%M%S)"
        mv "${LAI_RESULT}" "${backup}"
        echo "Existing LAI moved to: ${backup}"
    fi

    if [[ ! -s "${LAI_RESULT}" ]]; then
        echo "=== Calculate LAI: ${label} (${LAI_MODE}) ==="
        LAI -genome "${GENOME}" -intact "${INTACT}" -all "${ALL}" \
            -q -t "${THREADS}" "${LAI_EXTRA_ARGS[@]}"
    else
        echo "Existing LAI found; skipping: ${LAI_RESULT}"
    fi

    require_file "${LAI_RESULT}" "${label} LAI result"
    touch "${METHOD_MARKER}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' \
        "${label}" "${edta_dir}" "${lai_dir}" "${GENOME}" "${INTACT}" "${ALL}" \
        "${LAI_RESULT}" "${LAI_INPUT_METHOD}" "${LAI_MODE}" >> "${SUMMARY}"
done

echo
echo "=== All LAI calculations completed ==="
echo "Summary: ${SUMMARY}"
