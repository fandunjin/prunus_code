#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/karyotype}"
PREVIOUS_RIBBON_RUN="${PREVIOUS_RIBBON_RUN:-}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RIBBON_RUN_DIR:-${OUT_ROOT}/tree_pairwise_ribbon_${RUN_ID}}"
INPUT_DIR="${BASE_RUN}/inputs"
BASE_PAIR_DIR="${BASE_RUN}/pairs"
PAIR_DIR="${RUN_DIR}/anchors"
PREPARED_DIR="${RUN_DIR}/prepared"
FIG_DIR="${PREPARED_DIR}/figures"
GENOME_DIR="${GENOME_DIR:-${OUT_ROOT}/agr/genome}"
PREPARE_SCRIPT="${ROOT_DIR}/scripts/prepare_jcvi_tree_ribbon.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
PMU_MINSPAN="${PMU_MINSPAN:-30}"
MIN_CHR_MB="${MIN_CHR_MB:-5}"
MIN_BACKGROUND_SCORE="${MIN_BACKGROUND_SCORE:-0}"
MIN_CHR_PAIR_SCORE="${MIN_CHR_PAIR_SCORE:-0}"
MIN_CHR_PAIR_FRACTION="${MIN_CHR_PAIR_FRACTION:-0}"
MAX_BACKGROUND_BLOCKS_PER_CHR_PAIR="${MAX_BACKGROUND_BLOCKS_PER_CHR_PAIR:-0}"
HIGHLIGHT_ORDINALS="${HIGHLIGHT_ORDINALS:-12 15}"
HIGHLIGHT_COLORS="${HIGHLIGHT_COLORS:-#009E73 #E69F00}"
FIXED_ORDER_SPECIES="${FIXED_ORDER_SPECIES:-pcam}"
HIGHLIGHT_MAX_GENE_GAP="${HIGHLIGHT_MAX_GENE_GAP:-30}"
HIGHLIGHT_PAD_GENES="${HIGHLIGHT_PAD_GENES:-1}"
ORDER_SWEEPS="${ORDER_SWEEPS:-40}"
LOCAL_SWAP_PASSES="${LOCAL_SWAP_PASSES:-30}"
FIGSIZE="${FIGSIZE:-14x30}"

if [[ -z "$PREVIOUS_RIBBON_RUN" ]]; then
    shopt -s nullglob
    previous_candidates=("${OUT_ROOT}"/tree_pairwise_ribbon_*)
    shopt -u nullglob
    for candidate in "${previous_candidates[@]}"; do
        [[ -d "$candidate/anchors" ]] || continue
        if [[ -z "$PREVIOUS_RIBBON_RUN" || "$candidate" -nt "$PREVIOUS_RIBBON_RUN" ]]; then
            PREVIOUS_RIBBON_RUN="$candidate"
        fi
    done
fi
PREVIOUS_PAIR_DIR="${PREVIOUS_RIBBON_RUN:+${PREVIOUS_RIBBON_RUN}/anchors}"

# Leaf order read from the current 17-species tree.
TREE_ORDER=(pzi_copy1 psp ppa pmal pav pcam pspe ppe pdu psal pmu_refseq par msyl mdom pcom rch fve)

die() { echo "[ERROR] $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

[[ -s "$CONDA_SH" ]] || die "Conda initialization script not found: $CONDA_SH"
source "$CONDA_SH"
conda activate "$JCVI_ENV"
command -v python >/dev/null 2>&1 || die "python is unavailable in $JCVI_ENV"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in $JCVI_ENV"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in $JCVI_ENV"
python -c 'from ete4 import Tree' >/dev/null 2>&1 || die "ete4 is not importable in $JCVI_ENV"

for label in "${TREE_ORDER[@]}"; do
    require_file "$INPUT_DIR/${label}/${label}.bed"
    require_file "$INPUT_DIR/${label}/${label}.pep"
    require_file "$INPUT_DIR/${label}/${label}.seqids"
done

mkdir -p "$RUN_DIR" "$PAIR_DIR" "$FIG_DIR"
MANIFEST="$RUN_DIR/anchor_manifest.tsv"
: > "$MANIFEST"

echo "=== JCVI tree-ordered pairwise ribbon ==="
echo "Base run    : $BASE_RUN"
echo "Output      : $RUN_DIR"
echo "Reuse pairs : ${PREVIOUS_PAIR_DIR:-none}"
echo "Tree order  : ${TREE_ORDER[*]}"
echo "Min chr     : ${MIN_CHR_MB} Mb"
echo "Highlights  : pzi chromosomes ${HIGHLIGHT_ORDINALS} (${HIGHLIGHT_COLORS})"
echo "Fixed order : ${FIXED_ORDER_SPECIES}"
echo "Order sweeps: ${ORDER_SWEEPS}"
echo "Swap passes : ${LOCAL_SWAP_PASSES}"
echo "Background  : keep blocks with score >= ${MIN_BACKGROUND_SCORE}"
echo "Chr pair    : score >= ${MIN_CHR_PAIR_SCORE}; relative support >= ${MIN_CHR_PAIR_FRACTION}"
echo "PMU minspan : ${PMU_MINSPAN} (other pairs: ${MINSPAN})"
echo "Blocks/pair : ${MAX_BACKGROUND_BLOCKS_PER_CHR_PAIR} (0 means unlimited)"

for ((i=0; i<${#TREE_ORDER[@]}-1; i++)); do
    query="${TREE_ORDER[$i]}"
    ref="${TREE_ORDER[$((i+1))]}"
    comparison="${query}_vs_${ref}"
    work="$PAIR_DIR/$comparison"
    mkdir -p "$work"
    is_pmu_pair=0
    pair_minspan="$MINSPAN"
    if [[ "$query" == "pmu_refseq" || "$ref" == "pmu_refseq" ]]; then
        is_pmu_pair=1
        pair_minspan="$PMU_MINSPAN"
    fi

    source_simple="$BASE_PAIR_DIR/$comparison/${query}.${ref}.anchors.simple"
    if [[ "$is_pmu_pair" == "0" && -s "$source_simple" ]]; then
        echo "[INFO] Reuse existing anchors: $comparison"
        source_raw="$BASE_PAIR_DIR/$comparison/${query}.${ref}.anchors"
        [[ -s "$source_raw" ]] || source_raw="NA"
        printf '%s\t%s\t%s\t%s\n' "$query" "$ref" "$source_simple" "$source_raw" >> "$MANIFEST"
        continue
    fi

    previous_simple="${PREVIOUS_PAIR_DIR:-NA}/$comparison/${query}.${ref}.anchors.simple"
    previous_raw="${PREVIOUS_PAIR_DIR:-NA}/$comparison/${query}.${ref}.anchors"
    if [[ "$is_pmu_pair" == "0" && -n "${PREVIOUS_PAIR_DIR:-}" && -s "$previous_simple" ]]; then
        echo "[INFO] Reuse previous adjacent anchors: $comparison"
        [[ -s "$previous_raw" ]] || previous_raw="NA"
        printf '%s\t%s\t%s\t%s\n' "$query" "$ref" "$previous_simple" "$previous_raw" >> "$MANIFEST"
        continue
    fi

    for label in "$query" "$ref"; do
        ln -sfn "$INPUT_DIR/$label/$label.bed" "$work/$label.bed"
        ln -sfn "$INPUT_DIR/$label/$label.pep" "$work/$label.pep"
    done

    cd "$work"
    raw="${query}.${ref}.anchors"
    screened="${query}.${ref}.anchors.new"
    simple="${query}.${ref}.anchors.simple"
    if [[ "$is_pmu_pair" == "1" && -s "$previous_raw" && ! -s "$raw" ]]; then
        echo "[INFO] Re-screen existing PMU raw anchors with minspan=$pair_minspan: $comparison"
        ln -sfn "$previous_raw" "$raw"
    fi
    if [[ ! -s "$raw" ]]; then
        echo "[INFO] Catalog ortholog: $comparison"
        python -m jcvi.compara.catalog ortholog "$query" "$ref" \
            --dbtype=prot --align_soft=diamond_blastp \
            --cscore="$CSCORE" --min_size="$MIN_ANCHORS" \
            --cpus="$THREADS" --no_strip_names --no_dotplot
    fi
    require_file "$raw"
    if [[ ! -s "$screened" || ! -s "$simple" ]]; then
        python -m jcvi.compara.synteny screen \
            --qbed="$query.bed" --sbed="$ref.bed" \
            --minspan="$pair_minspan" --simple "$raw" "$screened"
    fi
    if [[ ! -s "$simple" && -s "$screened.simple" ]]; then
        ln -sfn "$screened.simple" "$simple"
    fi
    require_file "$simple"
    printf '%s\t%s\t%s\t%s\n' "$query" "$ref" "$work/$simple" "$work/$raw" >> "$MANIFEST"
    cd "$ROOT_DIR"
done

python "$PREPARE_SCRIPT" \
    --input-root "$INPUT_DIR" \
    --genome-dir "$GENOME_DIR" \
    --manifest "$MANIFEST" \
    --outdir "$PREPARED_DIR" \
    --species "${TREE_ORDER[@]}" \
    --pzi-label pzi_copy1 \
    --min-chr-mb "$MIN_CHR_MB" \
    --highlight-ordinals $HIGHLIGHT_ORDINALS \
    --highlight-colors $HIGHLIGHT_COLORS \
    --fixed-order-species $FIXED_ORDER_SPECIES \
    --highlight-max-gene-gap "$HIGHLIGHT_MAX_GENE_GAP" \
    --highlight-pad-genes "$HIGHLIGHT_PAD_GENES" \
    --order-sweeps "$ORDER_SWEEPS" \
    --local-swap-passes "$LOCAL_SWAP_PASSES" \
    --min-background-score "$MIN_BACKGROUND_SCORE" \
    --min-chr-pair-score "$MIN_CHR_PAIR_SCORE" \
    --min-chr-pair-fraction "$MIN_CHR_PAIR_FRACTION" \
    --max-background-blocks-per-chr-pair "$MAX_BACKGROUND_BLOCKS_PER_CHR_PAIR"

cd "$FIG_DIR"
python -m jcvi.graphics.karyotype \
    "$FIG_DIR/tree_order_pairwise.seqids" \
    "$FIG_DIR/tree_order_pairwise.layout" \
    --figsize="$FIGSIZE" \
    --format=png --notex --nocircles \
    --outfile="$FIG_DIR/tree_order_pairwise_ribbon.png"

python -m jcvi.graphics.karyotype \
    "$FIG_DIR/tree_order_pairwise.seqids" \
    "$FIG_DIR/tree_order_pairwise_basepair.layout" \
    --figsize="$FIGSIZE" \
    --format=png --notex --nocircles --basepair \
    --outfile="$FIG_DIR/tree_order_pairwise_ribbon_basepair.png"

echo "[OK] Tree-ordered pairwise ribbon: $FIG_DIR/tree_order_pairwise_ribbon.png"
echo "[OK] Base-pair-scaled ribbon: $FIG_DIR/tree_order_pairwise_ribbon_basepair.png"
echo "[OK] Base-pair track scale: $PREPARED_DIR/basepair_track_scale.tsv"
echo "[OK] Filter summary: $PREPARED_DIR/filtered_chromosomes.tsv"
echo "[OK] Highlight summary: $PREPARED_DIR/highlighted_pzi_chromosomes.tsv"
echo "[OK] Highlighted blocks: $PREPARED_DIR/highlighted_blocks.tsv"
echo "[OK] Order/orientation: $PREPARED_DIR/chromosome_order_orientation.tsv"
echo "[OK] Display labels: $PREPARED_DIR/chromosome_display_labels.tsv"
echo "[OK] Order metrics: $PREPARED_DIR/order_optimization.tsv"
echo "[OK] Background filter: $PREPARED_DIR/background_filtering.tsv"
