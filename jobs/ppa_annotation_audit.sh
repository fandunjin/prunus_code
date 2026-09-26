#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
CURRENT_POST="${CURRENT_POST:-${DATA_DIR}/annotation/eviann/postprocess}"
LEGACY_POST="${LEGACY_POST:-${DATA_DIR}/postprocess}"
EVIDENCE_DIR="${EVIDENCE_DIR:-${DATA_DIR}/annotation/eviann/evidence}"
HELIXER_GFF="${HELIXER_GFF:-${DATA_DIR}/annotation/helixer/ppa.helixer.gff3}"
EVIANN_RUN_DIR="${EVIANN_RUN_DIR:-${DATA_DIR}/annotation/eviann/run_ppa_unmasked}"
RAW_GFF="${RAW_GFF:-${EVIANN_RUN_DIR}/ppa.fa.pseudo_label.gff}"
GENOME="${GENOME:-${ROOT_DIR}/../prunus_old/analysis/prunus_pan_te/discovery/ppa/ppa.fa}"
OUT_DIR="${OUT_DIR:-${DATA_DIR}/annotation/eviann/audit}"
RUN_BUSCO="${RUN_BUSCO:-0}"
RUN_GENOME_BUSCO="${RUN_GENOME_BUSCO:-0}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
THREADS="${THREADS:-24}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }

count_fasta() {
    awk '/^>/{n++} END{print n+0}' "$1"
}

gff_counts() {
    awk -F '\t' '!/^#/ {if ($3=="gene") g++; else if ($3=="mRNA" || $3=="transcript") m++; else if ($3=="CDS") c++; else if ($3=="exon") e++} END {printf "%d\t%d\t%d\t%d\n", g+0,m+0,c+0,e+0}' "$1"
}

check_annotation() {
    local label="$1" post="$2" pep="${post}/ppa.eviann.protein.fa" gff="${post}/ppa.eviann.genes.gff3"
    if [[ ! -s "$pep" ]]; then
        printf '%s\tMISSING\t%s\t0\t0\t0\t0\t0\n' "$label" "$pep" >> "${OUT_DIR}/annotation_inventory.tsv"
        return 0
    fi
    local proteins genes mrnas cds exons
    proteins="$(count_fasta "$pep")"
    if [[ -s "$gff" ]]; then
        read -r genes mrnas cds exons <<< "$(gff_counts "$gff")"
    else
        genes=0; mrnas=0; cds=0; exons=0
    fi
    printf '%s\tOK\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$pep" "$proteins" "$genes" "$mrnas" "$cds" "$exons" \
        >> "${OUT_DIR}/annotation_inventory.tsv"
}

mkdir -p "$OUT_DIR"
printf 'label\tstatus\tprotein_fasta\tprotein_records\tgenes\tmRNAs\tCDS_parts\texons\n' \
    > "${OUT_DIR}/annotation_inventory.tsv"
check_annotation current "$CURRENT_POST"
check_annotation legacy "$LEGACY_POST"

{
    printf 'item\tpath\tstatus\n'
    for item_path in \
        "current_protein=${CURRENT_POST}/ppa.eviann.protein.fa" \
        "current_gff3=${CURRENT_POST}/ppa.eviann.genes.gff3" \
        "legacy_protein=${LEGACY_POST}/ppa.eviann.protein.fa" \
        "legacy_gff3=${LEGACY_POST}/ppa.eviann.genes.gff3" \
        "raw_eviann_gff=${RAW_GFF}" \
        "helixer_gff=${HELIXER_GFF}" \
        "eviann_unmasked_genome=${GENOME}" \
        "transcript_evidence=${EVIDENCE_DIR}/psp_pzi.cds.fa" \
        "protein_evidence=${EVIDENCE_DIR}/psp_pzi.proteins.fa"; do
        item="${item_path%%=*}"
        path="${item_path#*=}"
        [[ -s "$path" ]] && status=OK || status=MISSING
        printf '%s\t%s\t%s\n' "$item" "$path" "$status"
    done
} > "${OUT_DIR}/input_inventory.tsv"

if [[ -s "$RAW_GFF" ]]; then
    printf 'source\tgenes\tmRNAs\tCDS_parts\texons\n' > "${OUT_DIR}/gff_inventory.tsv"
    printf 'raw_eviann\t%s\n' "$(gff_counts "$RAW_GFF")" >> "${OUT_DIR}/gff_inventory.tsv"
fi
if [[ -s "$HELIXER_GFF" ]]; then
    [[ -s "${OUT_DIR}/gff_inventory.tsv" ]] || printf 'source\tgenes\tmRNAs\tCDS_parts\texons\n' > "${OUT_DIR}/gff_inventory.tsv"
    printf 'helixer\t%s\n' "$(gff_counts "$HELIXER_GFF")" >> "${OUT_DIR}/gff_inventory.tsv"
fi

for evidence in "${EVIDENCE_DIR}/psp_pzi.cds.fa" "${EVIDENCE_DIR}/psp_pzi.proteins.fa"; do
    if [[ -s "$evidence" ]]; then
        printf '%s\t%s\n' "$evidence" "$(count_fasta "$evidence")"
    fi
done > "${OUT_DIR}/evidence_record_counts.tsv"

if [[ "$RUN_BUSCO" == 1 ]]; then
    require_file "$CONDA_SH" "Conda initialization script"
    require_file "$BUSCO_LINEAGE" "BUSCO lineage dataset"
    source "$CONDA_SH"
    conda activate "$BUSCO_ENV"
    command -v busco >/dev/null 2>&1 || die "busco not found in ${BUSCO_ENV}"
    require_file "${CURRENT_POST}/ppa.eviann.protein.fa" "current PPA protein FASTA"
    busco -i "${CURRENT_POST}/ppa.eviann.protein.fa" -m proteins -l "$BUSCO_LINEAGE" \
        -c "$THREADS" --offline -o ppa_current_eviann_protein_eudicots \
        --out_path "$OUT_DIR/busco" --force
    if [[ -s "${LEGACY_POST}/ppa.eviann.protein.fa" ]]; then
        busco -i "${LEGACY_POST}/ppa.eviann.protein.fa" -m proteins -l "$BUSCO_LINEAGE" \
            -c "$THREADS" --offline -o ppa_legacy_eviann_protein_eudicots \
            --out_path "$OUT_DIR/busco" --force
    fi
fi

if [[ "$RUN_GENOME_BUSCO" == 1 ]]; then
    require_file "$CONDA_SH" "Conda initialization script"
    require_file "$BUSCO_LINEAGE" "BUSCO lineage dataset"
    require_file "$GENOME" "unmasked PPA genome"
    source "$CONDA_SH"
    conda activate "$BUSCO_ENV"
    command -v busco >/dev/null 2>&1 || die "busco not found in ${BUSCO_ENV}"
    busco -i "$GENOME" -m genome -l "$BUSCO_LINEAGE" -c "$THREADS" \
        --offline -o ppa_unmasked_genome_eudicots --out_path "$OUT_DIR/busco" --force
fi

cat > "${OUT_DIR}/README.txt" <<EOF
PPA annotation audit

annotation_inventory.tsv: current versus legacy protein and GFF3 counts.
input_inventory.tsv: existence of genome, Helixer, EviAnn and evidence files.
gff_inventory.tsv: raw EviAnn and Helixer feature counts.
evidence_record_counts.tsv: FASTA record counts for transcript/protein evidence.

The current postprocess script does not remove models; it copies the GFF3 and
exports sequences. A low current protein count therefore points upstream to
EviAnn input evidence/model generation, not to BUSCO settings.
EOF

echo "PPA audit completed: ${OUT_DIR}"
cat "${OUT_DIR}/annotation_inventory.tsv"
