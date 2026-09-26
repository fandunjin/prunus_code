#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# PZI copy1 self-synteny with the requested display order:
# 11,1,2,13,15,3,6,12,4,7,5,16,8,10,9,14.
# Only BED row order is changed; chromosome IDs, coordinates, genes and
# anchor relationships remain unchanged.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi}"
INPUT_DIR="${SOURCE_RUN}/inputs/pzi_copy1"
PAIR_DIR="${OUT_DIR}/pairs/pzi_copy1_self"
FIG_DIR="${OUT_DIR}/figures"
SUMMARY_DIR="${OUT_DIR}/summary"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"

LABEL="pzi_copy1"
SOURCE_BED="${INPUT_DIR}/${LABEL}.bed"
SOURCE_PEP="${INPUT_DIR}/${LABEL}.pep"
BED="${LABEL}.bed"
PEP="${LABEL}.pep"
RAW="${LABEL}.${LABEL}.anchors"
SCREENED="${LABEL}.${LABEL}.anchors.new"
SIMPLE="${LABEL}.${LABEL}.anchors.simple"

die() { echo "[ERROR] $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }
anchor_count() { awk 'NF >= 2 && $1 !~ /^#/ {n++} END {print n+0}' "$1"; }

mkdir -p "$PAIR_DIR" "$FIG_DIR" "$SUMMARY_DIR"
require_file "$CONDA_SH"; require_file "$SOURCE_BED"; require_file "$SOURCE_PEP"

# Reorder complete BED records by chromosome. Supports chr01 or C01-style IDs.
awk -F '\t' -v OFS='\t' '
    BEGIN {order="11 1 2 13 15 3 6 12 4 7 5 16 8 10 9 14"; n=split(order,a," ");
            for(i=1;i<=n;i++) rank[a[i]]=i }
    NF >= 4 {
        id=$1; gsub(/^C|^chr/, "", id); id=sprintf("%d", id+0)
        if (rank[id] > 0) {buf[rank[id], ++count[rank[id]]] = $0}
    }
    END {for(i=1;i<=n;i++) for(j=1;j<=count[i];j++) print buf[i,j]}
' "$SOURCE_BED" > "$PAIR_DIR/${BED}.tmp"
mv "$PAIR_DIR/${BED}.tmp" "$PAIR_DIR/${BED}"
ln -sfn "$SOURCE_PEP" "$PAIR_DIR/$PEP"

source "$CONDA_SH"
conda activate "$JCVI_ENV"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"
command -v python >/dev/null 2>&1 || die "python unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

if [[ "$FORCE" == "1" ]]; then
    rm -f "$PAIR_DIR/$RAW" "$PAIR_DIR/${RAW}.simple" "$PAIR_DIR/$SCREENED" \
        "$PAIR_DIR/${SCREENED}.simple" "$PAIR_DIR/$SIMPLE" \
        "$FIG_DIR/pzi_copy1_ordered_self.raw.dotplot.png" \
        "$FIG_DIR/pzi_copy1_ordered_self.raw.dotplot.pdf" \
        "$FIG_DIR/pzi_copy1_ordered_self.screened.dotplot.png" \
        "$FIG_DIR/pzi_copy1_ordered_self.screened.dotplot.pdf"
fi

cd "$PAIR_DIR"
if [[ ! -s "$RAW" ]]; then
    python -m jcvi.compara.catalog ortholog "$LABEL" "$LABEL" \
        --dbtype=prot --align_soft=diamond_blastp \
        --cscore="$CSCORE" --min_size="$MIN_ANCHORS" --cpus="$THREADS" \
        --no_strip_names --no_dotplot
fi
require_file "$RAW"
if [[ ! -s "$SCREENED" ]]; then
    python -m jcvi.compara.synteny screen --qbed="$BED" --sbed="$BED" \
        --minspan="$MINSPAN" --simple "$RAW" "$SCREENED"
fi
require_file "$SCREENED"
if [[ ! -s "$SIMPLE" ]]; then
    for f in "${SCREENED}.simple" "${RAW}.simple"; do
        if [[ -s "$f" ]]; then ln -sfn "$f" "$SIMPLE"; break; fi
    done
fi

for spec in \
    "$RAW|$FIG_DIR/pzi_copy1_ordered_self.raw.dotplot.png|png" \
    "$RAW|$FIG_DIR/pzi_copy1_ordered_self.raw.dotplot.pdf|pdf" \
    "$SCREENED|$FIG_DIR/pzi_copy1_ordered_self.screened.dotplot.png|png" \
    "$SCREENED|$FIG_DIR/pzi_copy1_ordered_self.screened.dotplot.pdf|pdf"; do
    IFS='|' read -r anchors output format <<< "$spec"
    if [[ ! -s "$output" ]]; then
        python -m jcvi.graphics.dotplot "$anchors" --qbed="$BED" --sbed="$BED" \
            --figsize=12x12 --dpi=300 --format="$format" --nochpf --notex --outfile="$output"
    fi
    require_file "$output"
done

printf 'label\torder\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\n' > "$SUMMARY_DIR/summary.tsv"
printf '%s\t11,1,2,13,15,3,6,12,4,7,5,16,8,10,9,14\t%s\t%s\t%s\t%s\n' \
    "$LABEL" "$(anchor_count "$RAW")" "$(anchor_count "$SCREENED")" \
    "$PAIR_DIR/$RAW" "$PAIR_DIR/$SCREENED" >> "$SUMMARY_DIR/summary.tsv"
echo "[OK] Ordered PZI copy1 self-dotplots completed: $FIG_DIR"
