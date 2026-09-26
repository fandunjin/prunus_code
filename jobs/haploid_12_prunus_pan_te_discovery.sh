#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# One independent EDTA discovery run per submitted target.  EDTA's structural
# discovery stages make effective use of this allocation; submit targets as
# separate jobs so the scheduler can place them independently.
# TARGET=all is supported for a sequential recovery run, but is not recommended.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/prunus_pan_te}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
THREADS="${THREADS:-24}"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
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
require_file "${EDTA_CDS}" "Shared peach CDS filter"
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v EDTA.pl >/dev/null 2>&1 || { echo "ERROR: EDTA.pl was not found." >&2; exit 1; }

if [ -z "${TARGET}" ]; then
    echo "ERROR: set TARGET to one of: ${PAN_TE_TARGETS[*]}" >&2
    exit 1
fi
if [ "${TARGET}" = "all" ]; then
    TARGETS=("${PAN_TE_TARGETS[@]}")
else
    TARGETS=("${TARGET}")
fi

mkdir -p "${OUT_ROOT}/discovery" "${LOG_DIR}"
echo "=== Standardized independent EDTA discovery ==="
echo "Targets     : ${TARGETS[*]}"
echo "Shared CDS  : ${EDTA_CDS}"
echo "EDTA env    : ${EDTA_ENV}"
echo "Threads     : ${THREADS}"
echo

for target in "${TARGETS[@]}"; do
    set_pan_te_target "${target}"
    require_file "${TARGET_SOURCE}" "Source genome for ${TARGET_LABEL}"
    work_dir="${OUT_ROOT}/discovery/${TARGET_LABEL}"
    work_fasta="${work_dir}/${TARGET_LABEL}.fa"
    library="${work_fasta}.mod.EDTA.TElib.fa"
    teanno="${work_fasta}.mod.EDTA.TEanno.gff3"
    summary="${work_fasta}.mod.EDTA.TEanno.sum"
    manifest="${work_dir}/${TARGET_LABEL}.discovery_manifest.tsv"
    mkdir -p "${work_dir}"

    if [ ! -e "${work_fasta}" ]; then
        if [ "${TARGET_GZIP}" = "1" ]; then
            echo ">>> Decompressing source genome: ${TARGET_LABEL}"
            gzip -dc "${TARGET_SOURCE}" > "${work_fasta}"
        else
            # panEDTA resolves symbolic links before looking for adjacent
            # *.mod.EDTA.* files.  Use a hard link where possible so the
            # staged basename remains stable; copy only across filesystems.
            echo ">>> Materializing source genome: ${TARGET_LABEL}"
            if ! ln "${TARGET_SOURCE}" "${work_fasta}" 2>/dev/null; then
                cp -p "${TARGET_SOURCE}" "${work_fasta}"
            fi
        fi
    fi
    require_file "${work_fasta}" "Working genome for ${TARGET_LABEL}"

    if [ ! -s "${teanno}" ]; then
        echo ">>> EDTA discovery: ${TARGET_LABEL}"
        # EDTA writes next to the process working directory rather than next
        # to an absolute --genome path.  Isolate every target so concurrent
        # jobs cannot put or mix their outputs in the project root.
        (
            cd "${work_dir}"
            EDTA.pl \
                --genome "$(basename "${work_fasta}")" \
                --species others \
                --step all \
                --anno 1 \
                --sensitive 1 \
                --force "${EDTA_FORCE}" \
                --overwrite "${EDTA_OVERWRITE}" \
                --threads "${THREADS}" \
                --cds "${EDTA_CDS}"
        )
    else
        echo ">>> Existing standardized EDTA result found: ${TARGET_LABEL}"
    fi
    require_file "${library}" "EDTA library for ${TARGET_LABEL}"
    require_file "${teanno}" "EDTA GFF3 for ${TARGET_LABEL}"
    require_file "${summary}" "EDTA summary for ${TARGET_LABEL}"
    {
        printf "target\tgroup\tsource_fasta\tworking_fasta\tshared_cds\tedta_library\tteanno_gff3\tteanno_summary\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TARGET_LABEL}" "${TARGET_GROUP}" "${TARGET_SOURCE}" "${work_fasta}" "${EDTA_CDS}" "${library}" "${teanno}" "${summary}"
    } > "${manifest}"
done

echo "=== Standardized EDTA discovery completed ==="
