#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/annotation/genes/pzi_eviann/postprocess/pzi32.eviann.genes.gff3}"
INPUT_ROOT="${INPUT_ROOT:-${ROOT_DIR}/haploid_repartition/inputs/pzi}"
OUT="${OUT:-${INPUT_ROOT}/pzi32_annotation_summary.tsv}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
count_gff() {
    local label="$1" gff="$2"
    require_file "$gff" "$label GFF3"
    awk -F '\t' -v label="$label" '
        BEGIN { gene=mrna=cds=exon=0 }
        !/^#/ && NF >= 9 {
            if ($3 == "gene") gene++
            else if ($3 == "mRNA" || $3 == "transcript") mrna++
            else if ($3 == "CDS") cds++
            else if ($3 == "exon") exon++
        }
        END { print label "\t" gene "\t" mrna "\t" cds "\t" exon }
    ' "$gff"
}

require_file "$COMPLETE_GFF" "PZI complete-assembly EviAnn GFF3"
mkdir -p "$(dirname "$OUT")" "${ROOT_DIR}/logs"
{
    printf 'Set\tGene number\tRNA number\tCDS number\tExon\n'
    count_gff pzi32 "$COMPLETE_GFF"
    for haplotype in hapA1 hapA2 hapB1 hapB2; do
        count_gff "$haplotype" "${INPUT_ROOT}/gff/pzi_${haplotype}.gff3"
    done
} > "${OUT}.tmp.$$"
mv -f "${OUT}.tmp.$$" "$OUT"
echo "[OK] PZI annotation summary: $OUT"
