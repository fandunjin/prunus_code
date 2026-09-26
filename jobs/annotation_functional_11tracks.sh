#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
FUNCTIONAL_ENV="${FUNCTIONAL_ENV:-orthofinder}"
DIAMOND_BIN="${DIAMOND_BIN:-diamond}"
SWISSPROT_DB="${SWISSPROT_DB:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
FROZEN_PEP="${FROZEN_PEP:-${ROOT_DIR}/annotation/functional/inputs_11tracks}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/functional}"
WORK_ROOT="${WORK_ROOT:-${OUT_ROOT}/work}"
THREADS="${THREADS:-24}"
SPECIES_SET="${SPECIES_SET:-pcam pav pspe pmal par pmu pdu psal ppe gtr rch}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || die "FORCE must be 0 or 1"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
require_file "${CONDA_SH}" "Conda initialization script"

source "${CONDA_SH}"
conda activate "${FUNCTIONAL_ENV}"
command -v "${DIAMOND_BIN}" >/dev/null 2>&1 || die "${DIAMOND_BIN} not found in ${FUNCTIONAL_ENV}"
"${DIAMOND_BIN}" version >/dev/null 2>&1 || die "${DIAMOND_BIN} failed to run in ${FUNCTIONAL_ENV}"

mkdir -p "${OUT_ROOT}" "${WORK_ROOT}" "${ROOT_DIR}/logs"

resolve_protein() {
    local track="$1"
    local candidate
    for candidate in "${FROZEN_PEP}/${track}.fa" \
                     "${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/pep/${track}.fa" \
                     "${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs/main/pep/${track}.fa" \
                     "${ROOT_DIR}/analysisphylogenetics/01_frozen_inputs/main/pep/${track}.fa"; do
        if [[ -s "${candidate}" ]]; then printf '%s\n' "${candidate}"; return 0; fi
    done
    die "protein FASTA for ${track} not found; stage the 11 species-level peps into ${FROZEN_PEP}"
}

prepare_diamond_db() {
    local source="$1"
    local label="$2"
    local db_prefix="${WORK_ROOT}/${label}"
    if [[ "$source" == *.dmnd ]]; then
        require_file "$source" "$label DIAMOND database"
        printf '%s\n' "$source"
        return 0
    fi
    require_file "$source" "$label database FASTA"
    if [[ ! -s "${db_prefix}.dmnd" ]]; then
        local input="$source"
        if [[ "$source" == *.gz ]]; then
            input="${WORK_ROOT}/${label}.fasta"
            if [[ ! -s "$input" ]]; then
                gzip -t "$source"
                gzip -dc "$source" > "${input}.tmp.$$"
                mv -f "${input}.tmp.$$" "$input"
            fi
        fi
        "${DIAMOND_BIN}" makedb --in "$input" --db "$db_prefix" \
            > "${db_prefix}.makedb.log" 2>&1
    fi
    require_file "${db_prefix}.dmnd" "$label DIAMOND database"
    printf '%s\n' "${db_prefix}.dmnd"
}

run_track() {
    local track="$1"
    local protein="$2"
    local out="${OUT_ROOT}/${track}"
    local result="${out}/swissprot.diamond.tsv"
    local normalized="${out}/protein_annotation.tsv"
    local go_table="${out}/go_annotation.tsv"
    local manifest="${out}/functional_annotation_manifest.txt"

    case "$track" in
        pzi|psp|ppa)
            die "refusing to annotate ${track}: racemose annotation is owned by jobs/annotation_functional_all.sh" ;;
    esac

    require_file "${protein}" "${track} frozen protein FASTA"
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
    require_file "${result}" "${track} Swiss-Prot DIAMOND result"

    {
        printf 'source\tqseqid\tsseqid\tpident\tlength\tevalue\tbitscore\tstitle\n'
        awk 'BEGIN {OFS="\t"} {print "swissprot",$1,$2,$3,$4,$11,$12,$13}' "${result}"
    } > "${normalized}"

    # Swiss-Prot DIAMOND output alone has no complete GO mapping; the metadata
    # transfer job fills swissprot_go.tsv / swissprot_functional_annotation.tsv.
    printf 'protein_id\tgo_terms\n' > "${go_table}"

    {
        echo "Species: ${track}"
        echo "Protein FASTA: ${protein}"
        echo "BLASTP/DIAMOND e-value: 1e-5"
        echo "BLASTP/DIAMOND max target sequences: 1"
        echo "Swiss-Prot enabled: 1 (${SWISSPROT_DB})"
        echo "NR enabled: 0"
        echo "KEGG enabled: 0"
        echo "KOG enabled: 0"
        echo "InterProScan enabled: 0"
        echo "DIAMOND environment: ${FUNCTIONAL_ENV}"
        echo "DIAMOND version: $(${DIAMOND_BIN} version 2>&1 | head -n 1)"
        echo "Outputs: ${normalized}; next step: jobs/annotation_functional_11tracks_metadata.sh"
        echo "Query FASTA scope: frozen OrthoFinder track (stable IDs join Orthogroups.tsv directly)"
    } > "${manifest}"

    echo "[OK] ${track}: swissprot hits=$(awk 'END {print NR-1}' "${normalized}")"
}

read -r -a TRACKS <<< "${SPECIES_SET//,/ }"
[[ "${#TRACKS[@]}" -gt 0 ]] || die "SPECIES_SET is empty"

for track in "${TRACKS[@]}"; do
    run_track "${track}" "$(resolve_protein "${track}")"
done

{
    echo "Functional annotation manifest (11 tracks)"
    echo "Tracks requested: ${SPECIES_SET}"
    echo "Queries: frozen OrthoFinder inputs under ${FROZEN_PEP}"
    echo "Database: ${SWISSPROT_DB}"
    echo "Method: DIAMOND blastp best hit, e-value 1e-5, max-target-seqs 1"
    echo "Outputs: ${OUT_ROOT}/<track>/swissprot.diamond.tsv + protein_annotation.tsv"
    echo "Next step: jsub < jobs/annotation_functional_11tracks_metadata.sh"
} > "${OUT_ROOT}/functional_annotation_11tracks_manifest.txt"

echo "[OK] All 11-track Swiss-Prot searches completed under ${OUT_ROOT}"
