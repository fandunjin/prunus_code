#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FUNCTIONAL_ENV="${FUNCTIONAL_ENV:-orthofinder}"
INTERPRO_ENV="${INTERPRO_ENV:-annot}"
THREADS="${THREADS:-24}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/functional}"
WORK_ROOT="${WORK_ROOT:-${OUT_ROOT}/work}"
DIAMOND_BIN="${DIAMOND_BIN:-diamond}"
INTERPROSCAN_BIN="${INTERPROSCAN_BIN:-interproscan.sh}"

# Swiss-Prot is the only database whose project path is known. Other databases
# are opt-in: set RUN_NR/RUN_KEGG/RUN_KOG to 1 and provide the corresponding
# FASTA or DIAMOND database path. Missing opt-in databases are fatal.
RUN_SWISSPROT="${RUN_SWISSPROT:-1}"
RUN_NR="${RUN_NR:-0}"
RUN_KEGG="${RUN_KEGG:-0}"
RUN_KOG="${RUN_KOG:-0}"
RUN_INTERPROSCAN="${RUN_INTERPROSCAN:-0}"
SWISSPROT_DB="${SWISSPROT_DB:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
# Place downloaded sequence databases under database/ with these names.  The
# corresponding RUN_* switches remain off by default until each local file is
# verified; the job accepts either FASTA(/.gz) or a prebuilt .dmnd database.
NR_DB="${NR_DB:-${ROOT_DIR}/database/NR/nr.gz}"
KEGG_DB="${KEGG_DB:-${ROOT_DIR}/database/KEGG/kegg.faa.gz}"
KOG_DB="${KOG_DB:-${ROOT_DIR}/database/KOG/kog.faa}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
check_bool() { [[ "$2" == 0 || "$2" == 1 ]] || die "$1 must be 0 or 1"; }

[[ -s "$CONDA_SH" ]] || die "Conda initialization script missing: $CONDA_SH"
check_bool RUN_SWISSPROT "$RUN_SWISSPROT"
check_bool RUN_NR "$RUN_NR"
check_bool RUN_KEGG "$RUN_KEGG"
check_bool RUN_KOG "$RUN_KOG"
check_bool RUN_INTERPROSCAN "$RUN_INTERPROSCAN"

source "$CONDA_SH"
conda activate "$FUNCTIONAL_ENV"
command -v "$DIAMOND_BIN" >/dev/null 2>&1 || die "diamond not found in ${FUNCTIONAL_ENV}"
if [[ "$RUN_INTERPROSCAN" == 1 ]]; then
    # InterProScan is maintained in the project's annotation environment,
    # whereas DIAMOND is maintained in the orthofinder environment.
    conda activate "$INTERPRO_ENV"
    command -v "$INTERPROSCAN_BIN" >/dev/null 2>&1 || die "InterProScan not found in ${INTERPRO_ENV}; set INTERPROSCAN_BIN"
    conda activate "$FUNCTIONAL_ENV"
fi
mkdir -p "$OUT_ROOT" "$WORK_ROOT" "${ROOT_DIR}/logs"

resolve_existing() {
    local value="$1"; shift
    if [[ -n "$value" ]]; then printf '%s' "$value"; return 0; fi
    local candidate
    for candidate in "$@"; do
        if [[ -s "$candidate" ]]; then printf '%s' "$candidate"; return 0; fi
    done
    printf '%s' "$1"
}

PZI_PROTEIN="${PZI_PROTEIN:-${ROOT_DIR}/annotation/genes/pzi_eviann/postprocess/pzi32.eviann.protein.fa}"
PSP_PROTEIN="${PSP_PROTEIN:-${ROOT_DIR}/data/psp/annotation/genes/psp32_eviann/postprocess/psp32.eviann.protein.fa}"
PPA_PROTEIN="${PPA_PROTEIN:-}"
PPA_PROTEIN="$(resolve_existing "$PPA_PROTEIN" \
  "${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.protein.fa" \
  "${ROOT_DIR}/data/ppa/postprocess/ppa.eviann.protein.fa" \
  "${ROOT_DIR}/data/ppa/annotation/genes/ppa_eviann/postprocess/ppa.eviann.protein.fa")"

prepare_diamond_db() {
    local source="$1"
    local label="$2"
    local db_prefix="${WORK_ROOT}/${label}"
    if [[ "$source" == *.dmnd ]]; then
        require_file "$source" "$label DIAMOND database"
        printf '%s' "$source"
        return 0
    fi
    require_file "$source" "$label database FASTA"
    if [[ ! -s "${db_prefix}.dmnd" ]]; then
        local input="$source"
        if [[ "$source" == *.gz ]]; then
            input="${WORK_ROOT}/${label}.fasta"
            if [[ ! -s "$input" ]]; then gzip -t "$source"; gzip -dc "$source" > "${input}.tmp.$$"; mv -f "${input}.tmp.$$" "$input"; fi
        fi
        "$DIAMOND_BIN" makedb --in "$input" --db "$db_prefix" > "${db_prefix}.makedb.log" 2>&1
    fi
    require_file "${db_prefix}.dmnd" "$label DIAMOND database"
    printf '%s' "${db_prefix}.dmnd"
}

run_diamond() {
    local species="$1"
    local protein="$2"
    local label="$3"
    local database="$4"
    local out="${OUT_ROOT}/${species}"
    local db_arg
    mkdir -p "$out"
    require_file "$protein" "${species} final protein FASTA"
    db_arg="$(prepare_diamond_db "$database" "$label")"
    local result="${out}/${label}.diamond.tsv"
    if [[ ! -s "$result" ]]; then
        "$DIAMOND_BIN" blastp --query "$protein" --db "$db_arg" --out "$result" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --evalue 1e-5 --max-target-seqs 1 --threads "$THREADS"
    fi
    require_file "$result" "${species} ${label} DIAMOND results"
    printf '%s\n' "$result"
}

run_species() {
    local species="$1"
    local protein="$2"
    local out="${OUT_ROOT}/${species}"
    local manifest="${out}/functional_annotation_manifest.txt"
    require_file "$protein" "${species} final protein FASTA"
    mkdir -p "$out"
    local result_files=()
    if [[ "$RUN_SWISSPROT" == 1 ]]; then result_files+=("$(run_diamond "$species" "$protein" swissprot "$SWISSPROT_DB")"); fi
    if [[ "$RUN_NR" == 1 ]]; then [[ -n "$NR_DB" ]] || die "RUN_NR=1 requires NR_DB"; result_files+=("$(run_diamond "$species" "$protein" nr "$NR_DB")"); fi
    if [[ "$RUN_KEGG" == 1 ]]; then [[ -n "$KEGG_DB" ]] || die "RUN_KEGG=1 requires KEGG_DB"; result_files+=("$(run_diamond "$species" "$protein" kegg "$KEGG_DB")"); fi
    if [[ "$RUN_KOG" == 1 ]]; then [[ -n "$KOG_DB" ]] || die "RUN_KOG=1 requires KOG_DB"; result_files+=("$(run_diamond "$species" "$protein" kog "$KOG_DB")"); fi

    if [[ "$RUN_INTERPROSCAN" == 1 ]]; then
        local base="${out}/interproscan"
        conda activate "$INTERPRO_ENV"
        if [[ ! -s "${base}.tsv" ]]; then
            "$INTERPROSCAN_BIN" -i "$protein" -b "$base" -f tsv,gff3 -goterms -pa -cpu "$THREADS" -iprlookup
        fi
        require_file "${base}.tsv" "${species} InterProScan TSV"
        [[ -s "${base}.gff3" ]] || die "${species} InterProScan GFF3 missing: ${base}.gff3"
        awk -F '\t' 'BEGIN{OFS="\t"; print "protein_id","go_terms"} $0 !~ /^#/ && NF >= 14 && $14 != "-" {print $1,$14}' "${base}.tsv" > "${out}/go_annotation.tsv"
        conda activate "$FUNCTIONAL_ENV"
    else
        printf 'protein_id\tgo_terms\n' > "${out}/go_annotation.tsv"
    fi

    {
        printf 'source\tqseqid\tsseqid\tpident\tlength\tevalue\tbitscore\tstitle\n'
        for result in "${result_files[@]}"; do
            source_name="$(basename "$result" .diamond.tsv)"
            awk -v source="$source_name" 'BEGIN{OFS="\t"} {print source,$1,$2,$3,$4,$11,$12,$13}' "$result"
        done
    } > "${out}/protein_annotation.tsv"
    {
        echo "Species: $species"
        echo "Protein FASTA: $protein"
        echo "BLASTP/DIAMOND e-value: 1e-5"
        echo "BLASTP/DIAMOND max target sequences: 1"
        echo "Swiss-Prot enabled: $RUN_SWISSPROT ($SWISSPROT_DB)"
        echo "NR enabled: $RUN_NR (${NR_DB:-not set})"
        echo "KEGG enabled: $RUN_KEGG (${KEGG_DB:-not set})"
        echo "KOG enabled: $RUN_KOG (${KOG_DB:-not set})"
        echo "InterProScan enabled: $RUN_INTERPROSCAN"
        echo "DIAMOND environment: $FUNCTIONAL_ENV"
        echo "InterProScan environment: $INTERPRO_ENV"
        echo "DIAMOND version: $($DIAMOND_BIN version 2>&1 | head -n 1)"
        echo "Outputs: ${out}/protein_annotation.tsv; ${out}/go_annotation.tsv"
        echo "Note: PPE/PCAM proteins used by EviAnn are structural evidence, not annotations of the target proteins."
    } > "$manifest"
    echo "${species} functional annotation completed: ${out}"
}

read -r -a SPECIES_LIST <<< "$SPECIES_SET"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty; choose one or more of: pzi psp ppa"
for species in "${SPECIES_LIST[@]}"; do
    case "$species" in
        pzi) run_species pzi "$PZI_PROTEIN" ;;
        psp) run_species psp "$PSP_PROTEIN" ;;
        ppa) run_species ppa "$PPA_PROTEIN" ;;
        *) die "Unknown species in SPECIES_SET: $species (choose pzi, psp or ppa)" ;;
    esac
done
{
    echo "Functional annotation manifest"
    echo "Species requested: $SPECIES_SET"
    echo "Databases are enabled explicitly through RUN_* variables."
    echo "InterProScan GO/domain output is disabled by default until its installation/database path is confirmed."
} > "${OUT_ROOT}/functional_annotation_manifest.txt"
echo "All functional annotations completed: ${OUT_ROOT}"
