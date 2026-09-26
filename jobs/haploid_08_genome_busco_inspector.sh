#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
GENOME_DIR="${HAP_RENAME_DIR}/genome"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
# Keep QC output locations stable; only the four-haplotype genome input moves
# to the renamed package.
QC_DIR="${ROOT_DIR}/haploid/final_haplotypes"
BUSCO_OUT="${BUSCO_OUT:-${QC_DIR}/busco_genome}"
INSPECTOR_OUT="${INSPECTOR_OUT:-${QC_DIR}/inspector}"
LOG_DIR="${LOG_DIR:-${QC_DIR}/logs}"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
INSPECTOR_ENV="${INSPECTOR_ENV:-inspector}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
ONT_READS="${ONT_READS:-${ROOT_DIR}/rawdata/HD-14-2DNA.pass.fastq.gz}"
HAP_READS_DIR="${HAP_READS_DIR:-}"
THREADS="${THREADS:-24}"
HAP_SET="${HAP_SET:-all}"
FORCE="${FORCE:-0}"
FORCE_INSPECTOR="${FORCE_INSPECTOR:-0}"

ALL_HAPS=(hapA1 hapA2 hapB1 hapB2)

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

case "${HAP_SET}" in
    all)
        HAPS=("${ALL_HAPS[@]}")
        ;;
    hapA1|hapA2|hapB1|hapB2)
        HAPS=("${HAP_SET}")
        ;;
    *)
        echo "ERROR: HAP_SET must be all, hapA1, hapA2, hapB1, or hapB2." >&2
        echo "Observed: ${HAP_SET}" >&2
        exit 1
        ;;
esac

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
require_file "${ONT_READS}" "ONT reads"
require_file "${PACKAGE_CHECKSUMS}" "renamed haplotype package checksums"
for hap in "${HAPS[@]}"; do
    require_file "${GENOME_DIR}/pzi_${hap}.fa" "${hap} genome FASTA"
done

mkdir -p "${BUSCO_OUT}" "${INSPECTOR_OUT}" "${LOG_DIR}"

set +u
source "${CONDA_SH}"
set -u

echo "=== PZI phased-haplotype genome QC ==="
echo "Haplotypes     : ${HAPS[*]}"
echo "Genome dir     : ${GENOME_DIR}"
echo "BUSCO output   : ${BUSCO_OUT}"
echo "Inspector out  : ${INSPECTOR_OUT}"
echo "Shared ONT     : ${ONT_READS}"
echo "Hap reads dir  : ${HAP_READS_DIR:-not set}"
echo "Threads        : ${THREADS}"
echo "Force BUSCO    : ${FORCE}"
echo "Force Inspector: ${FORCE_INSPECTOR}"
echo

echo "=== BUSCO genome mode ==="
conda activate "${BUSCO_ENV}"
command -v busco >/dev/null 2>&1 || {
    echo "ERROR: busco was not found in BUSCO_ENV=${BUSCO_ENV}." >&2
    exit 1
}

for hap in "${HAPS[@]}"; do
    genome="${GENOME_DIR}/pzi_${hap}.fa"
    run_name="pzi_${hap}_genome_eudicots"
    run_dir="${BUSCO_OUT}/${run_name}"
    saved_checksums="${run_dir}/.haploid_rename.checksums.md5"
    package_changed=0
    if ! cmp -s "${PACKAGE_CHECKSUMS}" "${saved_checksums}"; then
        package_changed=1
    fi
    existing_summary=$(find "${run_dir}" -maxdepth 3 -name 'short_summary*.txt' -type f -size +0c 2>/dev/null | head -n 1 || true)

    if [ -n "${existing_summary}" ] && [ "${FORCE}" != "1" ] && [ "${package_changed}" != "1" ]; then
        echo "Skip BUSCO ${hap}: existing summary ${existing_summary}"
        continue
    fi

    echo ">>> BUSCO ${hap}"
    busco_args=(
        -i "${genome}"
        -m genome
        -l "${BUSCO_LINEAGE}"
        -c "${THREADS}"
        --offline
        -o "${run_name}"
        --out_path "${BUSCO_OUT}"
    )
    if [ "${FORCE}" = "1" ] || [ "${package_changed}" = "1" ]; then
        busco_args+=(--force)
    elif [ -d "${run_dir}" ]; then
        echo "Restart incomplete BUSCO run: ${run_dir}"
        busco_args+=(--restart)
    fi
    busco "${busco_args[@]}"
    cp -f "${PACKAGE_CHECKSUMS}" "${saved_checksums}"
done

echo
echo "=== Inspector with ONT reads ==="
conda activate "${INSPECTOR_ENV}"
command -v inspector.py >/dev/null 2>&1 || {
    echo "ERROR: inspector.py was not found in INSPECTOR_ENV=${INSPECTOR_ENV}." >&2
    exit 1
}

for hap in "${HAPS[@]}"; do
    genome="${GENOME_DIR}/pzi_${hap}.fa"
    out_dir="${INSPECTOR_OUT}/${hap}"
    marker="${out_dir}/.inspector_completed"
    saved_checksums="${out_dir}/.haploid_rename.checksums.md5"
    reads="${ONT_READS}"
    package_changed=0
    if ! cmp -s "${PACKAGE_CHECKSUMS}" "${saved_checksums}"; then
        package_changed=1
    fi

    # Prefer pre-binned haplotype reads when available. Otherwise all ONT reads
    # are used, and cross-haplotype mappings must be considered during review.
    if [ -n "${HAP_READS_DIR}" ]; then
        for candidate in \
            "${HAP_READS_DIR}/pzi_${hap}.ont.fastq.gz" \
            "${HAP_READS_DIR}/pzi_${hap}.fastq.gz" \
            "${HAP_READS_DIR}/${hap}.fastq.gz"; do
            if [ -s "${candidate}" ]; then
                reads="${candidate}"
                break
            fi
        done
    fi
    require_file "${reads}" "Inspector reads for ${hap}"

    if [ -s "${marker}" ] && [ "${FORCE_INSPECTOR}" != "1" ] && [ "${package_changed}" != "1" ]; then
        echo "Skip Inspector ${hap}: completion marker exists: ${marker}"
        continue
    fi
    if [ -d "${out_dir}" ] && [ -n "$(find "${out_dir}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
        if [ "${FORCE_INSPECTOR}" = "1" ] || [ "${package_changed}" = "1" ]; then
            backup_dir="${out_dir}.legacy_$(date +%Y%m%d_%H%M%S)_$$"
            echo "[INFO] Preserve previous Inspector output: ${out_dir} -> ${backup_dir}"
            mv "${out_dir}" "${backup_dir}"
        else
            echo "ERROR: Inspector output is non-empty but has no completion marker: ${out_dir}" >&2
            echo "Move the incomplete directory aside or set FORCE_INSPECTOR=1 before rerunning ${hap}." >&2
            exit 1
        fi
    fi

    mkdir -p "${out_dir}"
    echo ">>> Inspector ${hap}"
    echo "Genome: ${genome}"
    echo "Reads : ${reads}"
    if [ "${reads}" = "${ONT_READS}" ]; then
        echo "WARNING: using the shared tetraploid ONT read set for ${hap}."
        echo "Cross-haplotype mappings can inflate Inspector error calls."
    fi

    inspector.py \
        -c "${genome}" \
        -r "${reads}" \
        -o "${out_dir}" \
        --datatype ont \
        -t "${THREADS}"

    printf 'completed\t%s\n' "$(date -Iseconds)" > "${marker}"
    cp -f "${PACKAGE_CHECKSUMS}" "${saved_checksums}"
done

echo
echo "=== Genome QC completed ==="
echo "BUSCO summaries:"
find "${BUSCO_OUT}" -maxdepth 4 -name 'short_summary*.txt' -type f -print | sort
echo
echo "Inspector completion markers:"
for hap in "${ALL_HAPS[@]}"; do
    marker="${INSPECTOR_OUT}/${hap}/.inspector_completed"
    if [ -s "${marker}" ]; then
        echo "${marker}"
    fi
done
