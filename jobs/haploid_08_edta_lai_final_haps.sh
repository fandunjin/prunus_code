#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# results below annotation/repeats/final_haplotypes for compatibility with the
# existing panEDTA/LAI jobs.
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
GENOME_DIR="${HAP_RENAME_DIR}/genome"
OUT_ROOT="${ROOT_DIR}/annotation/repeats/final_haplotypes"
LOG_DIR="${ROOT_DIR}/logs"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
EDTA_ENV="${EDTA_ENV:-EDTA}"
LAI_ENV="${LAI_ENV:-LAI}"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_lai_inputs.py"

EDTA_THREADS="${EDTA_THREADS:-24}"
LAI_THREADS="${LAI_THREADS:-24}"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
EDTA_MUTATION_RATE="${EDTA_MUTATION_RATE:-}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"

HAP_SET="${HAP_SET:-hapA1 hapA2 hapB1 hapB2}"
HAP_SET="${HAP_SET//,/ }"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

run_edta_for_hap() {
    local hap="$1"
    local source_fasta="$2"
    local edta_dir="$3"
    local work_fasta="${edta_dir}/pzi_${hap}.fa"
    local source_header work_header backup_dir

    mkdir -p "${edta_dir}"
    source_header="$(awk '/^>/{print; exit}' "${source_fasta}")"
    work_header=""
    if [ -s "${work_fasta}" ]; then
        work_header="$(awk '/^>/{print; exit}' "${work_fasta}")"
    fi
    # A previous run may have populated the compatibility output directory
    # with the old chr##_hap FASTA.  Do not attach those EDTA results to the
    # new C01--C08 package; move the stale run aside and start a fresh one.
    if [ -n "${work_header}" ] && [ "${source_header}" != "${work_header}" ]; then
        backup_dir="${edta_dir}.legacy_$(date +%Y%m%d_%H%M%S)_$$"
        echo "[INFO] Input chromosome IDs changed; moving stale EDTA run to ${backup_dir}"
        mv "${edta_dir}" "${backup_dir}"
        mkdir -p "${edta_dir}"
    fi
    if [ ! -s "${work_fasta}" ]; then
        cp -p "${source_fasta}" "${work_fasta}"
    fi

    source "${CONDA_SH}"
    conda activate "${EDTA_ENV}"
    cd "${edta_dir}"

    echo
    echo "=== EDTA repeat annotation: ${hap} ==="
    echo "Source FASTA : ${source_fasta}"
    echo "Work FASTA   : ${work_fasta}"
    echo "Out dir      : ${edta_dir}"
    echo "Threads      : ${EDTA_THREADS}"
    echo "Force        : ${EDTA_FORCE}"
    echo "Overwrite    : ${EDTA_OVERWRITE}"

    if [ -s "${EDTA_CDS}" ]; then
        echo "CDS filter   : ${EDTA_CDS}"
    else
        echo "CDS filter   : not used; CDS file was not found."
    fi

    if [ ! -s "${work_fasta}.mod.EDTA.TEanno.gff3" ]; then
        EDTA_CMD=(
            EDTA.pl
            --genome "${work_fasta}"
            --species others
            --step all
            --anno 1
            --sensitive 1
            --force "${EDTA_FORCE}"
            --overwrite "${EDTA_OVERWRITE}"
            --threads "${EDTA_THREADS}"
        )
        if [ -n "${EDTA_MUTATION_RATE}" ]; then
            echo "Mutation rate: ${EDTA_MUTATION_RATE}"
            EDTA_CMD+=(--u "${EDTA_MUTATION_RATE}")
        else
            echo "Mutation rate: not specified; EDTA will run without --u."
        fi
        if [ -s "${EDTA_CDS}" ]; then
            EDTA_CMD+=(--cds "${EDTA_CDS}")
        fi
        "${EDTA_CMD[@]}"
    else
        echo "Existing EDTA TE annotation found; skipping EDTA."
    fi

    require_file "${work_fasta}.mod.EDTA.TEanno.gff3" "${hap} EDTA TE annotation GFF3"
    require_file "${work_fasta}.mod.EDTA.TElib.fa" "${hap} EDTA TE library"
    require_file "${work_fasta}.mod.EDTA.TEanno.sum" "${hap} EDTA TE annotation summary"
}

run_lai_for_hap() {
    local hap="$1"
    local edta_dir="$2"
    local lai_dir="$3"
    local gff2bed_path="$4"
    local prefix="pzi_${hap}"

    mkdir -p "${lai_dir}"

    source "${CONDA_SH}"
    conda activate "${LAI_ENV}"

    if command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    else
        PYTHON_BIN="python"
    fi

    echo
    echo "=== Prepare LAI inputs: ${hap} ==="
    echo "EDTA dir : ${edta_dir}"
    echo "LAI dir  : ${lai_dir}"
    echo "gff2bed  : ${gff2bed_path}"

    "${PYTHON_BIN}" "${PREP_SCRIPT}" \
        --edta-dir "${edta_dir}" \
        --out-dir "${lai_dir}" \
        --prefix "${prefix}" \
        --gff2bed "${gff2bed_path}" \
        --method auto

    require_file "${lai_dir}/${prefix}.lai_inputs.env" "${hap} LAI env file"
    # shellcheck disable=SC1090
    source "${lai_dir}/${prefix}.lai_inputs.env"

    require_file "${GENOME}" "${hap} LAI genome"
    require_file "${INTACT}" "${hap} LAI intact LTR pass.list"
    require_file "${ALL}" "${hap} LAI all LTR RepeatMasker-like output"

    cd "${lai_dir}"

    echo
    echo "=== LAI calculation: ${hap} ==="
    echo "Genome : ${GENOME}"
    echo "Intact : ${INTACT}"
    echo "All    : ${ALL}"
    echo "Result : ${LAI_RESULT}"
    echo "Method : ${LAI_INPUT_METHOD}"
    echo "Threads: ${LAI_THREADS}"

    if [ ! -s "${LAI_RESULT}" ]; then
        LAI -genome "${GENOME}" -intact "${INTACT}" -all "${ALL}" -q -t "${LAI_THREADS}"
    else
        echo "Existing LAI output found; skipping LAI."
    fi

    require_file "${LAI_RESULT}" "${hap} LAI output"
}

mkdir -p "${OUT_ROOT}" "${LOG_DIR}"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${PREP_SCRIPT}" "LAI input preparation script"

if [ -n "${GFF2BED_SCRIPT}" ] && [ -s "${GFF2BED_SCRIPT}" ]; then
    GFF2BED_RESOLVED="${GFF2BED_SCRIPT}"
elif [ -s "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl" ]; then
    GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl"
elif [ -s "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl" ]; then
    GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl"
elif command -v gff2bed.pl >/dev/null 2>&1; then
    GFF2BED_RESOLVED="$(command -v gff2bed.pl)"
else
    echo "ERROR: gff2bed.pl was not found. Set GFF2BED_SCRIPT=/path/to/gff2bed.pl." >&2
    exit 1
fi

SUMMARY="${OUT_ROOT}/edta_lai_final_haplotypes_summary.tsv"
echo -e "hap_set\tfasta\tedta_dir\tteanno_gff3\tteanno_sum\tte_library\tlai_dir\tlai_output\tstatus" > "${SUMMARY}"

echo "=== EDTA + LAI for final haplotypes ==="
echo "Root      : ${ROOT_DIR}"
echo "Hap set   : ${HAP_SET}"
echo "Out root  : ${OUT_ROOT}"
echo "EDTA env  : ${EDTA_ENV}"
echo "LAI env   : ${LAI_ENV}"
echo "gff2bed   : ${GFF2BED_RESOLVED}"
echo

for hap in ${HAP_SET}; do
    case "${hap}" in
        hapA1|hapA2|hapB1|hapB2) ;;
        *)
            echo "ERROR: unsupported haplotype '${hap}'. Use hapA1, hapA2, hapB1, hapB2." >&2
            exit 1
            ;;
    esac

    HAP_FASTA="${GENOME_DIR}/pzi_${hap}.fa"
    if [ ! -s "${HAP_FASTA}" ]; then
        echo "ERROR: could not find FASTA for ${hap}" >&2
        echo "Checked: ${HAP_FASTA}" >&2
        exit 1
    fi

    EDTA_DIR="${OUT_ROOT}/${hap}_edta"
    LAI_DIR="${OUT_ROOT}/${hap}_lai"
    WORK_FASTA="${EDTA_DIR}/pzi_${hap}.fa"

    run_edta_for_hap "${hap}" "${HAP_FASTA}" "${EDTA_DIR}"
    run_lai_for_hap "${hap}" "${EDTA_DIR}" "${LAI_DIR}" "${GFF2BED_RESOLVED}"

    echo -e "${hap}\t${HAP_FASTA}\t${EDTA_DIR}\t${WORK_FASTA}.mod.EDTA.TEanno.gff3\t${WORK_FASTA}.mod.EDTA.TEanno.sum\t${WORK_FASTA}.mod.EDTA.TElib.fa\t${LAI_DIR}\t${LAI_DIR}/pzi_${hap}.all_ltr.out.LAI\tpass" >> "${SUMMARY}"

    echo
    echo "=== ${hap} EDTA + LAI completed ==="
done

echo
echo "=== All requested haplotypes completed ==="
echo "Summary: ${SUMMARY}"
