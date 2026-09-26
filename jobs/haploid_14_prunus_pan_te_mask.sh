#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/prunus_pan_te}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
PAN_VERSION="${PAN_VERSION:-prunus_pan_te_v1}"
PAN_LIBRARY="${PAN_LIBRARY:-${OUT_ROOT}/library/frozen/${PAN_VERSION}.fa}"
THREADS="${THREADS:-8}"
RUN_DIVERGENCE="${RUN_DIVERGENCE:-1}"
FORCE_RERUN="${FORCE_RERUN:-0}"
REPEATMASKER_CONFIG="${REPEATMASKER_CONFIG:-}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${TARGET_CONFIG}" "Pan-TE target configuration"
require_file "${PAN_LIBRARY}" "Frozen cross-species pan-TE library"
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v RepeatMasker >/dev/null 2>&1 || { echo "ERROR: RepeatMasker was not found." >&2; exit 1; }
command -v calcDivergenceFromAlign.pl >/dev/null 2>&1 || { echo "ERROR: calcDivergenceFromAlign.pl was not found." >&2; exit 1; }

# Conda commonly installs calcDivergenceFromAlign.pl in bin/ while keeping
# RepeatMaskerConfig.pm under share/RepeatMasker.  The utility does not always
# add that parent directory to @INC, so locate it and supply PERL5LIB explicitly.
if [ "${RUN_DIVERGENCE}" = "1" ]; then
    if [ -z "${REPEATMASKER_CONFIG}" ]; then
        for candidate in \
            "${CONDA_PREFIX}/share/RepeatMasker/RepeatMaskerConfig.pm" \
            "${CONDA_PREFIX}/lib/RepeatMasker/RepeatMaskerConfig.pm" \
            "$(dirname "$(readlink -f "$(command -v RepeatMasker)")")/RepeatMaskerConfig.pm"; do
            if [ -s "${candidate}" ]; then
                REPEATMASKER_CONFIG="${candidate}"
                break
            fi
        done
    fi
    if [ -z "${REPEATMASKER_CONFIG}" ]; then
        REPEATMASKER_CONFIG="$(find "${CONDA_PREFIX}" -type f -name RepeatMaskerConfig.pm -print -quit 2>/dev/null || true)"
    fi
    require_file "${REPEATMASKER_CONFIG}" "RepeatMaskerConfig.pm"
    REPEATMASKER_PERL5LIB="$(dirname "${REPEATMASKER_CONFIG}")"
    echo "RepeatMasker config: ${REPEATMASKER_CONFIG}"
fi

if [ -z "${TARGET}" ]; then
    echo "ERROR: set TARGET to one of: ${PAN_TE_TARGETS[*]}" >&2
    exit 1
fi
if [ "${TARGET}" = "all" ]; then
    TARGETS=("${PAN_TE_TARGETS[@]}")
else
    TARGETS=("${TARGET}")
fi

mkdir -p "${OUT_ROOT}/mask" "${LOG_DIR}"
for target in "${TARGETS[@]}"; do
    set_pan_te_target "${target}"
    require_file "${TARGET_SOURCE}" "Source genome for ${TARGET_LABEL}"
    mask_dir="${OUT_ROOT}/mask/${TARGET_LABEL}"
    input_fasta="${mask_dir}/${TARGET_LABEL}.fa"
    rmout="${input_fasta}.out"
    align="${input_fasta}.align"
    gff="${input_fasta}.out.gff"
    divergence="${input_fasta}.divsum"
    mkdir -p "${mask_dir}"
    if [ ! -e "${input_fasta}" ]; then
        if [ "${TARGET_GZIP}" = "1" ]; then
            gzip -dc "${TARGET_SOURCE}" > "${input_fasta}"
        else
            ln -s "${TARGET_SOURCE}" "${input_fasta}"
        fi
    fi
    require_file "${input_fasta}" "Masking input for ${TARGET_LABEL}"
    if [ "${FORCE_RERUN}" = "1" ]; then
        rm -f "${rmout}" "${align}" "${gff}" "${divergence}"
    fi
    if [ ! -s "${rmout}" ]; then
        echo ">>> RepeatMasker standardized annotation: ${TARGET_LABEL}"
        RepeatMasker \
            -pa "${THREADS}" \
            -lib "${PAN_LIBRARY}" \
            -s \
            -nolow \
            -norna \
            -no_is \
            -a \
            -gff \
            -dir "${mask_dir}" \
            "${input_fasta}"
    else
        echo ">>> Existing standardized RepeatMasker output found: ${TARGET_LABEL}"
    fi
    require_file "${rmout}" "RepeatMasker .out for ${TARGET_LABEL}"
    require_file "${align}" "RepeatMasker .align for ${TARGET_LABEL}"
    if [ "${RUN_DIVERGENCE}" = "1" ] && [ ! -s "${divergence}" ]; then
        echo ">>> RepeatMasker divergence calculation: ${TARGET_LABEL}"
        PERL5LIB="${REPEATMASKER_PERL5LIB}${PERL5LIB:+:${PERL5LIB}}" \
            calcDivergenceFromAlign.pl -s "${divergence}" "${align}"
    fi
    printf "target\tgroup\tpan_library\tmask_input_fasta\trepeatmasker_out\trepeatmasker_align\tdivergence_summary\n" > "${mask_dir}/${TARGET_LABEL}.mask_manifest.tsv"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${TARGET_LABEL}" "${TARGET_GROUP}" "${PAN_LIBRARY}" "${input_fasta}" "${rmout}" "${align}" "${divergence}" \
        >> "${mask_dir}/${TARGET_LABEL}.mask_manifest.tsv"
done

echo "=== Standardized pan-library masking completed ==="
