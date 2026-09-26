#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Re-run Swiss-Prot functional searches on each already extracted haplotype
# protein FASTA.  Structural annotation and gene prediction are not rerun.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/repartition/functional}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
FUNCTIONAL_ENV="${FUNCTIONAL_ENV:-orthofinder}"
DIAMOND_BIN="${DIAMOND_BIN:-diamond}"
SWISSPROT_DB="${SWISSPROT_DB:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
THREADS="${THREADS:-24}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SWISSPROT_DB}" "Swiss-Prot FASTA or DIAMOND database"
[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || die "FORCE must be 0 or 1"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"

source "${CONDA_SH}"
conda activate "${FUNCTIONAL_ENV}"
command -v "${DIAMOND_BIN}" >/dev/null 2>&1 || die "${DIAMOND_BIN} not found in ${FUNCTIONAL_ENV}"
mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"

haplotypes_for() {
    case "$1" in
        pzi|psp) printf '%s\n' hapA1 hapA2 hapB1 hapB2 ;;
        ppa) printf '%s\n' hapA hapB ;;
        *) die "Unknown species: $1 (choose pzi, psp or ppa)" ;;
    esac
}

prepare_diamond_db() {
    local source="$1"
    local label="$2"
    local db_prefix="${OUT_ROOT}/work/${label}"
    mkdir -p "${OUT_ROOT}/work"
    if [[ "${source}" == *.dmnd ]]; then
        printf '%s\n' "${source}"
        return 0
    fi
    if [[ ! -s "${db_prefix}.dmnd" ]]; then
        local input="${source}"
        if [[ "${source}" == *.gz ]]; then
            input="${db_prefix}.fasta"
            if [[ ! -s "${input}" ]]; then
                gzip -t "${source}"
                gzip -dc "${source}" > "${input}.tmp.$$"
                mv -f "${input}.tmp.$$" "${input}"
            fi
        fi
        "${DIAMOND_BIN}" makedb --in "${input}" --db "${db_prefix}" \
            > "${db_prefix}.makedb.log" 2>&1
    fi
    require_file "${db_prefix}.dmnd" "Swiss-Prot DIAMOND database"
    printf '%s\n' "${db_prefix}.dmnd"
}

run_one() {
    local species="$1" hap="$2"
    local prefix="${species}_${hap}"
    local protein="${INPUT_ROOT}/${species}/pep/${prefix}.pep.fa"
    local out="${OUT_ROOT}/${species}/${hap}"
    local result="${out}/${prefix}.swissprot.diamond.tsv"
    local normalized="${out}/${prefix}.protein_annotation.tsv"
    local go_table="${out}/${prefix}.go_annotation.tsv"
    local manifest="${out}/${prefix}.functional_annotation_manifest.txt"
    require_file "${protein}" "${species}/${hap} protein FASTA"
    mkdir -p "${out}"
    if [[ "${FORCE}" == 1 ]]; then
        rm -f "${result}" "${normalized}" "${go_table}" "${manifest}"
    fi
    local db
    db="$(prepare_diamond_db "${SWISSPROT_DB}" swissprot)"
    if [[ ! -s "${result}" ]]; then
        "${DIAMOND_BIN}" blastp --query "${protein}" --db "${db}" --out "${result}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
            sstart send evalue bitscore stitle --evalue 1e-5 --max-target-seqs 1 \
            --threads "${THREADS}"
    fi
    require_file "${result}" "${species}/${hap} Swiss-Prot DIAMOND result"
    {
        printf 'source\tqseqid\tsseqid\tpident\tlength\tevalue\tbitscore\tstitle\n'
        awk 'BEGIN {OFS="\t"} {print "swissprot",$1,$2,$3,$4,$11,$12,$13}' "${result}"
    } > "${normalized}"
    # Swiss-Prot DIAMOND output does not provide a complete GO mapping. Keep a
    # stable empty table so downstream summaries distinguish this from missing.
    printf 'protein_id\tgo_terms\n' > "${go_table}"
    {
        echo "Species: ${species}"
        echo "Haplotype: ${hap}"
        echo "Protein FASTA: ${protein}"
        echo "Database: ${SWISSPROT_DB}"
        echo "DIAMOND e-value: 1e-5"
        echo "DIAMOND max target sequences: 1"
        echo "Threads: ${THREADS}"
        echo "Scope: haplotype-specific independent Swiss-Prot search"
        echo "Output: ${normalized}"
        echo "GO annotation: not inferred from Swiss-Prot-only search"
    } > "${manifest}"
    echo "[OK] ${species}/${hap}: Swiss-Prot hits=$(awk 'END {print NR-1}' "${normalized}")"
}

read -r -a SPECIES_LIST <<< "${SPECIES_SET//,/ }"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"
for species in "${SPECIES_LIST[@]}"; do
    while IFS= read -r hap; do
        run_one "${species}" "${hap}"
    done < <(haplotypes_for "${species}")
done

echo "[OK] Haplotype-specific Swiss-Prot annotation completed under ${OUT_ROOT}"
