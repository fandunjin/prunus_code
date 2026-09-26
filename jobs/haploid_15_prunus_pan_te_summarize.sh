#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
PAN_VERSION="${PAN_VERSION:-prunus_pan_te_v1}"
PAN_LIBRARY="${PAN_LIBRARY:-${OUT_ROOT}/library/frozen/${PAN_VERSION}.fa}"
THREADS="${THREADS:-8}"
FORCE_RERUN="${FORCE_RERUN:-0}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SUMMARY_SCRIPT="${SUMMARY_SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_summary.py}"

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
require_file "${SUMMARY_SCRIPT}" "Pan-TE summary helper"
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 was not found." >&2; exit 1; }
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

summary_dir="${OUT_ROOT}/summary"
mapping_dir="${summary_dir}/original_to_pan_family"
targets_tsv="${summary_dir}/standardized_mask_targets.tsv"
mkdir -p "${summary_dir}" "${mapping_dir}" "${LOG_DIR}"
printf "target\tgroup\tmask_input_fasta\trepeatmasker_out\tdiscovery_library\n" > "${targets_tsv}"

for target in "${PAN_TE_TARGETS[@]}"; do
    set_pan_te_target "${target}"
    discovery_fasta="${OUT_ROOT}/discovery/${TARGET_LABEL}/${TARGET_LABEL}.fa"
    discovery_library="${discovery_fasta}.mod.EDTA.TElib.fa"
    mask_dir="${OUT_ROOT}/mask/${TARGET_LABEL}"
    mask_input="${mask_dir}/${TARGET_LABEL}.fa"
    rmout="${mask_input}.out"
    paf="${mapping_dir}/${TARGET_LABEL}.discovery_to_${PAN_VERSION}.paf"
    map_tsv="${mapping_dir}/${TARGET_LABEL}.original_to_pan_family.tsv"
    require_file "${discovery_library}" "Discovery library for ${TARGET_LABEL}"
    require_file "${mask_input}" "Standardized masking input for ${TARGET_LABEL}"
    require_file "${rmout}" "Standardized RepeatMasker output for ${TARGET_LABEL}"
    if [ "${FORCE_RERUN}" = "1" ]; then
        rm -f "${paf}" "${map_tsv}"
    fi
    if [ ! -s "${paf}" ]; then
        echo ">>> Mapping original EDTA library to frozen pan library: ${TARGET_LABEL}"
        minimap2 -x asm20 --secondary=no -t "${THREADS}" "${PAN_LIBRARY}" "${discovery_library}" > "${paf}"
    fi
    if [ ! -s "${map_tsv}" ]; then
        "${PYTHON_BIN}" "${SUMMARY_SCRIPT}" map-library \
            --label "${TARGET_LABEL}" \
            --original-lib "${discovery_library}" \
            --paf "${paf}" \
            --output "${map_tsv}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "${TARGET_LABEL}" "${TARGET_GROUP}" "${mask_input}" "${rmout}" "${discovery_library}" \
        >> "${targets_tsv}"
done

echo ">>> Summarizing unified RepeatMasker annotations"
"${PYTHON_BIN}" "${SUMMARY_SCRIPT}" summarize \
    --targets-tsv "${targets_tsv}" \
    --out-dir "${summary_dir}"

printf "Pan-TE summary outputs\n======================\n\n" > "${summary_dir}/review_checklist.txt"
printf -- "- Use prunus_pan_te_total_summary.tsv for total masked-repeat content.\n" >> "${summary_dir}/review_checklist.txt"
printf -- "- Use prunus_pan_te_class_summary.tsv and prunus_pan_te_family_summary.tsv for the cross-species main comparisons.\n" >> "${summary_dir}/review_checklist.txt"
printf -- "- Use original_to_pan_family/*.tsv to audit classification migration from the standardized independent EDTA libraries.\n" >> "${summary_dir}/review_checklist.txt"
printf -- "- Preserve the independent EDTA intact-element outputs for structural TE and insertion-age analyses.\n" >> "${summary_dir}/review_checklist.txt"

echo "=== Cross-species pan-TE summary completed ==="
