#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Calculate LAI from the already repartitioned EDTA outputs. EDTA is not rerun.
# Native .real.pass.list/.all_ltr.out files are preferred; otherwise
# prepare_lai_inputs.py rebuilds them from the EDTA TEanno GFF3.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
LAI_ROOT="${LAI_ROOT:-${ROOT_DIR}/annotation/repartition/lai}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
LAI_THREADS="${LAI_THREADS:-24}"
LAI_METHOD="${LAI_METHOD:-auto}"
FORCE_PREP="${FORCE_PREP:-0}"
FORCE_LAI="${FORCE_LAI:-0}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
LAI_ENV="${LAI_ENV:-LAI}"
PREP_SCRIPT="${PREP_SCRIPT:-${ROOT_DIR}/scripts/prepare_lai_inputs.py}"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() {
    local path="$1" label="$2"
    [[ -s "${path}" ]] || die "Missing or empty ${label}: ${path}"
}

resolve_gff2bed() {
    if [[ -n "${GFF2BED_SCRIPT}" && -s "${GFF2BED_SCRIPT}" ]]; then
        printf '%s\n' "${GFF2BED_SCRIPT}"
    elif [[ -s "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl" ]]; then
        printf '%s\n' "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl"
    elif [[ -s "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl" ]]; then
        printf '%s\n' "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl"
    elif command -v gff2bed.pl >/dev/null 2>&1; then
        command -v gff2bed.pl
    else
        die "gff2bed.pl was not found; set GFF2BED_SCRIPT=/path/to/gff2bed.pl"
    fi
}

haplotypes_for() {
    case "$1" in
        pzi|psp) printf '%s\n' hapA1 hapA2 hapB1 hapB2 ;;
        ppa) printf '%s\n' hapA hapB ;;
        *) die "Unknown species: $1 (choose pzi, psp or ppa)" ;;
    esac
}

extract_score() {
    awk 'BEGIN { FS="\t" } $1 == "whole_genome" { print $7; found=1; exit } END { if (!found) print "NA" }' "$1"
}

extract_genome_bp() {
    awk 'BEGIN { FS="\t" } NR == 2 { print $3; found=1; exit } END { if (!found) print "NA" }' "$1"
}

run_one() {
    local species="$1" hap="$2"
    local edta_dir="${EDTA_INPUT_ROOT}/${species}/edta/${hap}"
    local prefix="${species}_${hap}"
    local lai_dir="${LAI_ROOT}/${species}/${hap}_lai"
    local te_gff3="${edta_dir}/${prefix}.fa.mod.EDTA.TEanno.gff3"
    local te_sum="${edta_dir}/${prefix}.fa.mod.EDTA.TEanno.sum"
    local env_file="${lai_dir}/${prefix}.lai_inputs.env"
    local genome intact all result backup score genome_bp method

    # LAI changes the working directory; restore the project root before each
    # target so user-supplied relative paths remain stable across targets.
    cd "${ROOT_DIR}"
    require_file "${te_gff3}" "${species}/${hap} EDTA TE annotation GFF3"
    require_file "${te_sum}" "${species}/${hap} EDTA TE summary"
    mkdir -p "${lai_dir}"

    if [[ "${FORCE_PREP}" == "1" ]]; then
        rm -f "${env_file}" "${lai_dir}/${prefix}.lai_inputs.tsv" \
            "${lai_dir}/${prefix}.fa.EDTA.TEanno.struc.bed" \
            "${lai_dir}/${prefix}.EDTA.TEanno.struc.bed"
    fi
    if [[ "${FORCE_LAI}" == "1" ]]; then
        # LAI creates these intermediates beside the final .LAI table.  Remove
        # stale files so a rerun cannot reuse an empty LTR FASTA/LTRlist.
        rm -f "${lai_dir}/${prefix}.all_ltr.out.LAI" \
            "${lai_dir}/${prefix}.all_ltr.out.LAI."*
    fi

    echo
    echo "=== Prepare LAI inputs: ${species}/${hap} ==="
    echo "EDTA dir : ${edta_dir}"
    echo "LAI dir  : ${lai_dir}"
    "${PYTHON_BIN}" "${PREP_SCRIPT}" \
        --edta-dir "${edta_dir}" \
        --out-dir "${lai_dir}" \
        --prefix "${prefix}" \
        --gff2bed "${GFF2BED_RESOLVED}" \
        --method "${LAI_METHOD}"

    require_file "${env_file}" "${species}/${hap} LAI input environment"
    # shellcheck disable=SC1090
    source "${env_file}"
    genome="${GENOME}"
    intact="${INTACT}"
    all="${ALL}"
    result="${LAI_RESULT}"
    method="${LAI_INPUT_METHOD}"
    require_file "${genome}" "${species}/${hap} LAI genome"
    require_file "${intact}" "${species}/${hap} intact LTR input"
    require_file "${all}" "${species}/${hap} all-LTR input"

    if [[ "${FORCE_LAI}" == "1" && -s "${result}" ]]; then
        backup="${result}.previous_$(date +%Y%m%d_%H%M%S)_$$"
        mv "${result}" "${backup}"
        echo "[INFO] Existing LAI result moved to: ${backup}"
    fi

    cd "${lai_dir}"
    echo "=== Calculate LAI: ${species}/${hap} ==="
    echo "Genome : ${genome}"
    echo "Intact : ${intact}"
    echo "All    : ${all}"
    echo "Method : ${method}"
    echo "Result : ${result}"
    echo "Threads: ${LAI_THREADS}"
    if [[ ! -s "${result}" ]]; then
        LAI -genome "${genome}" -intact "${intact}" -all "${all}" \
            -q -t "${LAI_THREADS}"
    else
        echo "[SKIP] Existing LAI result: ${result}"
    fi
    require_file "${result}" "${species}/${hap} LAI result"
    # LAI can still write a numeric .LAI table after failing to extract the
    # LTR sequences.  Do not accept that silent partial result.
    if [[ -e "${lai_dir}/${prefix}.all_ltr.out.LAI.LTR.fa" &&
          ! -s "${lai_dir}/${prefix}.all_ltr.out.LAI.LTR.fa" ]]; then
        die "${species}/${hap}: LAI produced an empty LTR FASTA; check genome/LAI IDs"
    fi

    score="$(extract_score "${result}")"
    genome_bp="$(extract_genome_bp "${result}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' \
        "${species}" "${hap}" "${edta_dir}" "${lai_dir}" "${genome}" \
        "${intact}" "${all}" "${result}" "${method}" "${genome_bp}" \
        "${score}" >> "${SUMMARY_TMP}"
    printf '%s\t%s\n' "${species}/${hap}" "${score}" >> "${SCORES_TMP}"
    echo "[OK] ${species}/${hap} LAI: ${score}"
}

require_file "${CONDA_SH}" "conda initialization script"
require_file "${PREP_SCRIPT}" "LAI input preparation script"
[[ "${FORCE_PREP}" == "0" || "${FORCE_PREP}" == "1" ]] || die "FORCE_PREP must be 0 or 1"
[[ "${FORCE_LAI}" == "0" || "${FORCE_LAI}" == "1" ]] || die "FORCE_LAI must be 0 or 1"
[[ "${LAI_THREADS}" =~ ^[1-9][0-9]*$ ]] || die "LAI_THREADS must be a positive integer"

source "${CONDA_SH}"
conda activate "${LAI_ENV}"
command -v LAI >/dev/null 2>&1 || die "LAI is unavailable in ${LAI_ENV}"
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable in ${LAI_ENV}"
PYTHON_VERSION="$("${PYTHON_BIN}" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
    || die "Unable to query Python version: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' \
    || die "Python >= 3.6 is required; found ${PYTHON_VERSION}"
GFF2BED_RESOLVED="$(resolve_gff2bed)"
mkdir -p "${LAI_ROOT}" "${ROOT_DIR}/logs"

SPECIES_SET="${SPECIES_SET//,/ }"
read -r -a SPECIES_LIST <<< "${SPECIES_SET}"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"

SUMMARY="${LAI_ROOT}/haplotype_lai_summary.tsv"
SUMMARY_TMP="${SUMMARY}.tmp.$$"
SCORES_TMP="${LAI_ROOT}/haplotype_lai_scores.tsv.tmp.$$"
trap 'rm -f "${SUMMARY_TMP}" "${SCORES_TMP}"' EXIT
printf 'species\thaplotype\tedta_dir\tlai_dir\tgenome\tintact_ltr\tall_ltr\tlai_result\tinput_method\tgenome_bp\twhole_genome_lai\tstatus\n' > "${SUMMARY_TMP}"
printf 'species_haplotype\twhole_genome_lai\n' > "${SCORES_TMP}"

echo "=== Haplotype LAI calculation ==="
echo "Species   : ${SPECIES_SET}"
echo "EDTA root : ${EDTA_INPUT_ROOT}"
echo "LAI root  : ${LAI_ROOT}"
echo "LAI env   : ${LAI_ENV}"
echo "Python    : $(command -v "${PYTHON_BIN}") (${PYTHON_VERSION})"
echo "Threads   : ${LAI_THREADS}"
echo "Method    : ${LAI_METHOD}"
echo "Force prep: ${FORCE_PREP}"
echo "Force LAI : ${FORCE_LAI}"
echo "gff2bed   : ${GFF2BED_RESOLVED}"

for species in "${SPECIES_LIST[@]}"; do
    while IFS= read -r hap; do
        run_one "${species}" "${hap}"
    done < <(haplotypes_for "${species}")
done

{
    head -n 1 "${SUMMARY_TMP}"
    tail -n +2 "${SUMMARY_TMP}" | sort -t $'\t' -k1,1 -k2,2
} > "${SUMMARY}.sorted"
mv "${SUMMARY}.sorted" "${SUMMARY}"
mv "${SCORES_TMP}" "${LAI_ROOT}/haplotype_lai_scores.tsv"
echo
echo "[OK] LAI summary: ${SUMMARY}"
echo "[OK] LAI scores : ${LAI_ROOT}/haplotype_lai_scores.tsv"
