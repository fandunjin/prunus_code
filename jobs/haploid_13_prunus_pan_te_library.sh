#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
PAN_VERSION="${PAN_VERSION:-prunus_pan_te_v1}"
# Retain single-copy full-length candidates during pan-library construction.
# A threshold of 3, appropriate for a PZI-only within-individual library, can
# discard true TE families restricted to one Prunus lineage or haplotype.
FULL_LENGTH_COPY="${FULL_LENGTH_COPY:-1}"
THREADS="${THREADS:-24}"
FORCE_RERUN="${FORCE_RERUN:-0}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PANEDTA_BIN="${PANEDTA_BIN:-panEDTA.sh}"

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
command -v "${PANEDTA_BIN}" >/dev/null 2>&1 || { echo "ERROR: ${PANEDTA_BIN} was not found." >&2; exit 1; }

library_dir="${OUT_ROOT}/library"
draft_dir="${library_dir}/draft"
frozen_dir="${library_dir}/frozen"
genome_list="${library_dir}/prunus_pan_te_genomes.list"
manifest="${library_dir}/prunus_pan_te_discovery_manifest.tsv"
draft_library="${draft_dir}/${PAN_VERSION}.draft.fa"
frozen_library="${frozen_dir}/${PAN_VERSION}.fa"
mkdir -p "${library_dir}" "${draft_dir}" "${frozen_dir}" "${LOG_DIR}"

: > "${genome_list}"
printf "target\tgroup\tworking_fasta\tedta_library\n" > "${manifest}"
for target in "${PAN_TE_TARGETS[@]}"; do
    set_pan_te_target "${target}"
    work_fasta="${OUT_ROOT}/discovery/${TARGET_LABEL}/${TARGET_LABEL}.fa"
    library="${work_fasta}.mod.EDTA.TElib.fa"
    require_file "${work_fasta}" "Standardized discovery genome ${TARGET_LABEL}"
    if [ -L "${work_fasta}" ]; then
        echo "ERROR: ${work_fasta} is a symbolic link." >&2
        echo "panEDTA resolves symlinks and may rerun EDTA in the source directory." >&2
        echo "Run scripts/materialize_prunus_pan_te_fastas.sh first." >&2
        exit 1
    fi
    require_file "${library}" "Standardized discovery library ${TARGET_LABEL}"
    printf "%s\n" "${work_fasta}" >> "${genome_list}"
    printf "%s\t%s\t%s\t%s\n" "${TARGET_LABEL}" "${TARGET_GROUP}" "${work_fasta}" "${library}" >> "${manifest}"
done

if [ "${FORCE_RERUN}" = "1" ]; then
    rm -f "${draft_library}" "${frozen_library}" "${frozen_library}.sha256"
fi

if [ ! -s "${draft_library}" ]; then
    echo "=== Building ${PAN_VERSION} draft pan-TE library ==="
    cd "${library_dir}"
    "${PANEDTA_BIN}" \
        -g "${genome_list}" \
        -c "${EDTA_CDS}" \
        -t "${THREADS}" \
        -f "${FULL_LENGTH_COPY}" \
        -a 0 \
        -o 0
    generated_library="${genome_list}.panEDTA.TElib.fa"
    require_file "${generated_library}" "panEDTA draft library"
    cp -p "${generated_library}" "${draft_library}"
else
    echo "=== Existing draft pan-TE library found ==="
fi

require_file "${draft_library}" "Draft pan-TE library"
if command -v TEsorter >/dev/null 2>&1; then
    if [ ! -s "${draft_dir}/${PAN_VERSION}.draft.fa.rexdb-plant.cls.tsv" ]; then
        echo ">>> TEsorter audit of draft pan-TE library"
        TEsorter "${draft_library}" -db rexdb-plant -p "${THREADS}"
    fi
else
    echo "WARNING: TEsorter is unavailable; classification audit is skipped." >&2
fi

# This copy freezes the exact machine-generated library used for standardized
# masking. Inspect the TEsorter audit before treating a questionable family as
# a biological conclusion; do not overwrite this version after masking starts.
if [ ! -s "${frozen_library}" ]; then
    cp -p "${draft_library}" "${frozen_library}"
    sha256sum "${frozen_library}" > "${frozen_library}.sha256"
fi

require_file "${frozen_library}" "Frozen pan-TE library"
printf "library_version\tlibrary\tsha256_file\tfull_length_copy\tshared_cds\ttarget_count\n" > "${library_dir}/frozen_library_manifest.tsv"
printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${PAN_VERSION}" "${frozen_library}" "${frozen_library}.sha256" "${FULL_LENGTH_COPY}" "${EDTA_CDS}" "${#PAN_TE_TARGETS[@]}" \
    >> "${library_dir}/frozen_library_manifest.tsv"

echo "=== Cross-species pan-TE library completed ==="
echo "Frozen library: ${frozen_library}"
