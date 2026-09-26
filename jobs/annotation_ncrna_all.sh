#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
TRNASCAN_ENV="${TRNASCAN_ENV:-infernal}"
INFERNAL_ENV="${INFERNAL_ENV:-infernal}"
THREADS="${THREADS:-16}"
# The Rfam database is prepared once with `cmpress -F` and kept under the
# project so batch jobs do not depend on an interactive shell variable.
RFAM_CM="${RFAM_CM:-${ROOT_DIR}/database/Rfam/Rfam.cm}"
# Avoid carrying the documentation placeholder from an earlier interactive
# command into a batch submission.
[[ "$RFAM_CM" == "/path/to/Rfam.cm" ]] && RFAM_CM="${ROOT_DIR}/database/Rfam/Rfam.cm"
TRNASCAN_BIN="${TRNASCAN_BIN:-tRNAscan-SE}"
CMSCAN_BIN="${CMSCAN_BIN:-cmscan}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/ncrna}"
WORK_ROOT="${WORK_ROOT:-${ROOT_DIR}/annotation/ncrna/work}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }

[[ -s "$CONDA_SH" ]] || die "Conda initialization script missing: $CONDA_SH"
require_file "$RFAM_CM" "Rfam covariance-model database"
for suffix in i1f i1i i1m i1p; do
    require_file "${RFAM_CM}.${suffix}" "Rfam pressed index (${suffix}); run cmpress -F ${RFAM_CM} in the infernal environment"
done

source "$CONDA_SH"
conda activate "$TRNASCAN_ENV"
command -v "$TRNASCAN_BIN" >/dev/null 2>&1 || die "tRNAscan-SE not found in ${TRNASCAN_ENV}"
conda activate "$INFERNAL_ENV"
command -v "$CMSCAN_BIN" >/dev/null 2>&1 || die "cmscan not found in ${INFERNAL_ENV}"
mkdir -p "$OUT_ROOT" "$WORK_ROOT" "${ROOT_DIR}/logs"

prepare_genome() {
    local source="$1" target="$2"
    if [[ -s "$target" ]]; then return 0; fi
    mkdir -p "$(dirname "$target")"
    if [[ "$source" == *.gz ]]; then
        require_file "$source" "compressed genome source"
        gzip -t "$source"
        gzip -dc "$source" > "${target}.tmp.$$"
    else
        require_file "$source" "genome source"
        cp -p "$source" "${target}.tmp.$$"
    fi
    mv -f "${target}.tmp.$$" "$target"
}

run_species() {
    local species="$1" source_genome="$2"
    local genome="${WORK_ROOT}/${species}.genome.fa"
    local out="${OUT_ROOT}/${species}"
    local trna_out="${out}/tRNAscan-SE.out"
    local trna_gff="${out}/tRNAscan-SE.gff3"
    local rfam_tbl="${out}/rfam.cmscan.tblout"
    local rfam_gff="${out}/rfam.cmscan.gff3"
    local rfam_log="${out}/rfam.cmscan.log"
    prepare_genome "$source_genome" "$genome"
    mkdir -p "$out"

    if [[ ! -s "$trna_out" || ! -s "$trna_gff" ]]; then
        conda activate "$TRNASCAN_ENV"
        "$TRNASCAN_BIN" --thread "$THREADS" --output "$trna_out" --gff "$trna_gff" "$genome"
    else
        echo "Existing tRNAscan output found for ${species}; skipping"
    fi
    require_file "$trna_out" "${species} tRNAscan output"
    require_file "$trna_gff" "${species} tRNAscan GFF3"
    awk '!/^#/ && NF == 9 {ok=1} END {exit(ok ? 0 : 1)}' "$trna_gff" \
        || die "${species} tRNAscan output is not valid 9-column GFF3: ${trna_gff}"

    if [[ ! -s "$rfam_tbl" || ! -s "$rfam_gff" ]]; then
        conda activate "$INFERNAL_ENV"
        # Infernal format 2 is the standard GFF3 main output; --tblout is
        # retained separately for reproducible hit filtering and auditing.
        "$CMSCAN_BIN" --cpu "$THREADS" --rfam --cut_ga --fmt 2 --tblout "$rfam_tbl" \
            "$RFAM_CM" "$genome" > "$rfam_gff" 2> "$rfam_log"
    else
        echo "Existing Infernal/Rfam output found for ${species}; skipping"
    fi
    require_file "$rfam_tbl" "${species} Rfam cmscan tblout"
    require_file "$rfam_gff" "${species} Rfam cmscan GFF3"
    awk '!/^#/ && NF == 9 {ok=1} END {exit(ok ? 0 : 1)}' "$rfam_gff" \
        || die "${species} cmscan output is not valid 9-column GFF3; check Infernal --fmt support"
    combined_gff="${out}/ncrna.combined.gff3"
    if [[ ! -s "$combined_gff" ]]; then
        {
            echo '##gff-version 3'
            grep -v '^#' "$trna_gff"
            grep -v '^#' "$rfam_gff"
        } > "${combined_gff}.tmp.$$"
        mv -f "${combined_gff}.tmp.$$" "$combined_gff"
    fi
    require_file "$combined_gff" "${species} combined ncRNA GFF3"

    trna_count="$(awk '!/^#/ && NF == 9 {n++} END {print n+0}' "$trna_gff")"
    rfam_count="$(awk '!/^#/ && NF == 9 {n++} END {print n+0}' "$rfam_gff")"
    {
        printf 'species\tclass\tcount\tmethod\n'
        printf '%s\ttRNA\t%s\ttRNAscan-SE v2\n' "$species" "$trna_count"
        printf '%s\tRfam_cm_hit\t%s\tInfernal cmscan + Rfam\n' "$species" "$rfam_count"
    } > "${out}/ncrna_summary.tsv"
    {
        echo "Species: $species"
        echo "Genome: $genome"
        echo "Rfam CM: $RFAM_CM"
        echo "tRNAscan-SE: $TRNASCAN_BIN"
        echo "Infernal cmscan: $CMSCAN_BIN"
        echo "tRNAscan-SE environment: $TRNASCAN_ENV"
        echo "Infernal environment: $INFERNAL_ENV"
        echo "Threads: $THREADS"
        echo "tRNAscan-SE version: $($TRNASCAN_BIN --version 2>&1 | head -n 1)"
        echo "Infernal version: $($CMSCAN_BIN -h 2>&1 | head -n 1)"
        echo "tRNA output: $trna_out"
        echo "Rfam tblout: $rfam_tbl"
        echo "Rfam GFF3: $rfam_gff"
        echo "Combined ncRNA GFF3: $combined_gff"
        echo "Note: Rfam hits are reported as covariance-model matches and require duplicate/overlap review before a final non-coding RNA catalogue is released."
    } > "${out}/ncrna_manifest.txt"
    echo "${species} non-coding RNA annotation completed: ${out}"
}

PZI_GENOME="${PZI_GENOME:-${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa}"
PSP_GENOME="${PSP_GENOME:-${ROOT_DIR}/data/psp/haphic/merge/04.build/final_assembly_qc_outputs/psp_2n32.chromosomes.fa}"
PPA_GENOME="${PPA_GENOME:-${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
read -r -a SPECIES_LIST <<< "$SPECIES_SET"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty; choose one or more of pzi psp ppa"
for species in "${SPECIES_LIST[@]}"; do
    case "$species" in
        pzi) run_species pzi "$PZI_GENOME" ;;
        psp) run_species psp "$PSP_GENOME" ;;
        ppa) run_species ppa "$PPA_GENOME" ;;
        *) die "Unknown species in SPECIES_SET: $species (choose pzi, psp or ppa)" ;;
    esac
done

{
    echo "Non-coding RNA annotation manifest"
    echo "Rfam CM: $RFAM_CM"
    echo "tRNAscan-SE environment: $TRNASCAN_ENV"
    echo "Infernal environment: $INFERNAL_ENV"
    echo "Species requested: ${SPECIES_SET}"
    for species in "${SPECIES_LIST[@]}"; do echo "${species}: ${OUT_ROOT}/${species}"; done
    echo "tRNA method: tRNAscan-SE v2"
    echo "rRNA/snRNA/miRNA method: Infernal cmscan against Rfam covariance models"
} > "${OUT_ROOT}/ncrna_manifest.txt"
echo "All non-coding RNA annotations completed: ${OUT_ROOT}"
