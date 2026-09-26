#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Re-run ncRNA discovery on each already extracted haplotype genome.  This is
# intentionally independent of the complete-assembly ncRNA job; no EDTA or
# gene annotation program is run here.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/repartition/ncrna}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
TRNASCAN_ENV="${TRNASCAN_ENV:-infernal}"
INFERNAL_ENV="${INFERNAL_ENV:-infernal}"
TRNASCAN_BIN="${TRNASCAN_BIN:-tRNAscan-SE}"
CMSCAN_BIN="${CMSCAN_BIN:-cmscan}"
RFAM_CM="${RFAM_CM:-${ROOT_DIR}/database/Rfam/Rfam.cm}"
THREADS="${THREADS:-16}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${RFAM_CM}" "Rfam covariance-model database"
for suffix in i1f i1i i1m i1p; do
    require_file "${RFAM_CM}.${suffix}" "Rfam pressed index ${suffix}"
done
[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || die "FORCE must be 0 or 1"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"

source "${CONDA_SH}"
conda activate "${TRNASCAN_ENV}"
command -v "${TRNASCAN_BIN}" >/dev/null 2>&1 || die "${TRNASCAN_BIN} not found in ${TRNASCAN_ENV}"
conda activate "${INFERNAL_ENV}"
command -v "${CMSCAN_BIN}" >/dev/null 2>&1 || die "${CMSCAN_BIN} not found in ${INFERNAL_ENV}"
mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"

haplotypes_for() {
    case "$1" in
        pzi|psp) printf '%s\n' hapA1 hapA2 hapB1 hapB2 ;;
        ppa) printf '%s\n' hapA hapB ;;
        *) die "Unknown species: $1 (choose pzi, psp or ppa)" ;;
    esac
}

run_one() {
    local species="$1" hap="$2"
    local prefix="${species}_${hap}"
    local genome="${INPUT_ROOT}/${species}/genome/${prefix}.fa"
    local out="${OUT_ROOT}/${species}/${hap}"
    local trna_out="${out}/${prefix}.tRNAscan-SE.out"
    local trna_gff="${out}/${prefix}.tRNAscan-SE.gff3"
    local rfam_tbl="${out}/${prefix}.rfam.cmscan.tblout"
    local rfam_gff="${out}/${prefix}.rfam.cmscan.gff3"
    local rfam_log="${out}/${prefix}.rfam.cmscan.log"
    local combined="${out}/${prefix}.ncrna.combined.gff3"

    require_file "${genome}" "${species}/${hap} genome"
    mkdir -p "${out}"
    if [[ "${FORCE}" == 1 ]]; then
        rm -f "${trna_out}" "${trna_gff}" "${rfam_tbl}" "${rfam_gff}" \
            "${rfam_log}" "${combined}" "${out}/${prefix}.ncrna_summary.tsv" \
            "${out}/${prefix}.ncrna_manifest.txt"
    fi

    if [[ ! -s "${trna_out}" || ! -s "${trna_gff}" ]]; then
        conda activate "${TRNASCAN_ENV}"
        "${TRNASCAN_BIN}" --thread "${THREADS}" --output "${trna_out}" \
            --gff "${trna_gff}" "${genome}"
    fi
    require_file "${trna_out}" "${species}/${hap} tRNAscan output"
    require_file "${trna_gff}" "${species}/${hap} tRNAscan GFF3"
    awk '!/^#/ && NF == 9 {ok=1} END {exit(ok ? 0 : 1)}' "${trna_gff}" \
        || die "Invalid tRNAscan GFF3: ${trna_gff}"

    if [[ ! -s "${rfam_tbl}" || ! -s "${rfam_gff}" ]]; then
        conda activate "${INFERNAL_ENV}"
        "${CMSCAN_BIN}" --cpu "${THREADS}" --rfam --cut_ga --fmt 2 \
            --tblout "${rfam_tbl}" "${RFAM_CM}" "${genome}" \
            >"${rfam_gff}" 2>"${rfam_log}"
    fi
    require_file "${rfam_tbl}" "${species}/${hap} Rfam tblout"
    require_file "${rfam_gff}" "${species}/${hap} Rfam GFF3"
    awk '!/^#/ && NF == 9 {ok=1} END {exit(ok ? 0 : 1)}' "${rfam_gff}" \
        || die "Invalid cmscan GFF3: ${rfam_gff}"

    if [[ ! -s "${combined}" ]]; then
        {
            echo '##gff-version 3'
            echo "# Independent ncRNA scan of ${prefix} haplotype genome"
            grep -v '^#' "${trna_gff}"
            grep -v '^#' "${rfam_gff}"
        } > "${combined}.tmp.$$"
        mv -f "${combined}.tmp.$$" "${combined}"
    fi
    require_file "${combined}" "${species}/${hap} combined ncRNA GFF3"

    local trna_count rfam_count
    trna_count="$(awk '!/^#/ && NF == 9 {n++} END {print n+0}' "${trna_gff}")"
    rfam_count="$(awk '!/^#/ && NF == 9 {n++} END {print n+0}' "${rfam_gff}")"
    {
        printf 'species\thaplotype\tclass\tcount\tmethod\n'
        printf '%s\t%s\ttRNA\t%s\ttRNAscan-SE v2\n' "${species}" "${hap}" "${trna_count}"
        printf '%s\t%s\tRfam_cm_hit\t%s\tInfernal cmscan + Rfam\n' "${species}" "${hap}" "${rfam_count}"
    } > "${out}/${prefix}.ncrna_summary.tsv"
    {
        echo "Species: ${species}"
        echo "Haplotype: ${hap}"
        echo "Genome: ${genome}"
        echo "Rfam CM: ${RFAM_CM}"
        echo "Threads: ${THREADS}"
        echo "tRNA GFF3: ${trna_gff}"
        echo "Rfam GFF3: ${rfam_gff}"
        echo "Combined ncRNA GFF3: ${combined}"
        echo "Scope: haplotype-specific independent scan"
    } > "${out}/${prefix}.ncrna_manifest.txt"
    echo "[OK] ${species}/${hap}: tRNA=${trna_count}, Rfam=${rfam_count}"
}

read -r -a SPECIES_LIST <<< "${SPECIES_SET//,/ }"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"
for species in "${SPECIES_LIST[@]}"; do
    while IFS= read -r hap; do
        run_one "${species}" "${hap}"
    done < <(haplotypes_for "${species}")
done

echo "[OK] Haplotype-specific ncRNA annotation completed under ${OUT_ROOT}"
