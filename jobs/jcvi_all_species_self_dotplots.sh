#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Run the same JCVI self-synteny workflow and parameters for every species
# already prepared by the main JCVI run. Each species keeps its existing BED
# chromosome names and protein IDs.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
SOURCE_INPUT_ROOT="${SOURCE_INPUT_ROOT:-${SOURCE_RUN}/inputs}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/all_species_self_jcvi}"
PAIR_ROOT="${OUT_DIR}/pairs"
FIG_DIR="${OUT_DIR}/figures"
SUMMARY_DIR="${OUT_DIR}/summary"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
# Deliberately identical to the PCAM/PPE/GTR script.
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }
anchor_count() { awk 'NF >= 2 && $1 !~ /^#/ {n++} END {print n+0}' "$1"; }

mkdir -p "$PAIR_ROOT" "$FIG_DIR" "$SUMMARY_DIR"
require_file "$CONDA_SH"
[[ -d "$SOURCE_INPUT_ROOT" ]] || die "Input directory not found: $SOURCE_INPUT_ROOT"

source "$CONDA_SH"
conda activate "$JCVI_ENV"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

printf 'species\tchromosomes\tgenes\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot_png\tscreened_dotplot_png\n' \
    > "${SUMMARY_DIR}/self_synteny_summary.tsv"

run_self() {
    local species="$1"
    local input_dir="${SOURCE_INPUT_ROOT}/${species}"
    local pair_dir="${PAIR_ROOT}/${species}_self"
    local bed="${species}.bed"
    local pep="${species}.pep"
    local raw="${species}.${species}.anchors"
    local screened="${species}.${species}.anchors.new"
    local simple="${species}.${species}.anchors.simple"
    local raw_png="${FIG_DIR}/${species}_self.raw.dotplot.png"
    local raw_pdf="${FIG_DIR}/${species}_self.raw.dotplot.pdf"
    local screened_png="${FIG_DIR}/${species}_self.screened.dotplot.png"
    local screened_pdf="${FIG_DIR}/${species}_self.screened.dotplot.pdf"

    require_file "${input_dir}/${bed}"
    require_file "${input_dir}/${pep}"
    mkdir -p "$pair_dir"
    ln -sfn "${input_dir}/${bed}" "${pair_dir}/${bed}"
    ln -sfn "${input_dir}/${pep}" "${pair_dir}/${pep}"

    if [[ "$FORCE" == "1" ]]; then
        rm -f "${pair_dir}/${raw}" "${pair_dir}/${raw}.simple" \
            "${pair_dir}/${screened}" "${pair_dir}/${screened}.simple" \
            "${pair_dir}/${simple}" "$raw_png" "$raw_pdf" "$screened_png" "$screened_pdf"
    fi

    (
        cd "$pair_dir"
        if [[ ! -s "$raw" ]]; then
            echo "[RUN] ${species}: catalog ortholog self"
            python -m jcvi.compara.catalog ortholog "$species" "$species" \
                --dbtype=prot --align_soft=diamond_blastp \
                --cscore="$CSCORE" --min_size="$MIN_ANCHORS" \
                --cpus="$THREADS" --no_strip_names --no_dotplot
        else
            echo "[SKIP] ${species}: existing raw anchors"
        fi
        require_file "$raw"

        if [[ ! -s "$screened" ]]; then
            python -m jcvi.compara.synteny screen \
                --qbed="$bed" --sbed="$bed" --minspan="$MINSPAN" \
                --simple "$raw" "$screened"
        else
            echo "[SKIP] ${species}: existing screened anchors"
        fi
        require_file "$screened"

        if [[ ! -s "$simple" ]]; then
            for candidate in "${screened}.simple" "${raw}.simple"; do
                if [[ -s "$candidate" ]]; then ln -sfn "$candidate" "$simple"; break; fi
            done
        fi

        for spec in "${raw}|${raw_png}|png" "${raw}|${raw_pdf}|pdf" \
                    "${screened}|${screened_png}|png" "${screened}|${screened_pdf}|pdf"; do
            IFS='|' read -r anchors output format <<< "$spec"
            if [[ ! -s "$output" ]]; then
                python -m jcvi.graphics.dotplot "$anchors" \
                    --qbed="$bed" --sbed="$bed" --figsize=12x12 --dpi=300 \
                    --format="$format" --nochpf --notex --outfile="$output"
            else
                echo "[SKIP] ${species}: existing figure $output"
            fi
            require_file "$output"
        done
    )

    local nchr ngenes
    nchr="$(awk 'NF >= 4 {seen[$1]=1} END {for (x in seen) n++; print n+0}' "${input_dir}/${bed}")"
    ngenes="$(awk 'NF >= 4 {n++} END {print n+0}' "${input_dir}/${bed}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species" "$nchr" "$ngenes" "$(anchor_count "${pair_dir}/${raw}")" \
        "$(anchor_count "${pair_dir}/${screened}")" "${pair_dir}/${raw}" \
        "${pair_dir}/${screened}" "$raw_png" "$screened_png" \
        >> "${SUMMARY_DIR}/self_synteny_summary.tsv"
}

# Use every prepared species directory that contains both canonical JCVI files.
mapfile -t SPECIES < <(
    find "$SOURCE_INPUT_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
        | while read -r s; do
            [[ -s "${SOURCE_INPUT_ROOT}/${s}/${s}.bed" && -s "${SOURCE_INPUT_ROOT}/${s}/${s}.pep" ]] && echo "$s"
          done | sort
)
[[ "${#SPECIES[@]}" -gt 0 ]] || die "No prepared species with BED and PEP files found"

echo "=== All-species JCVI self-synteny ==="
echo "Input root : ${SOURCE_INPUT_ROOT}"
echo "Output     : ${OUT_DIR}"
echo "Species    : ${SPECIES[*]}"
echo "JCVI env   : ${JCVI_ENV}"
echo "DIAMOND    : $(command -v diamond)"
echo "Parameters : cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}"

for species in "${SPECIES[@]}"; do
    run_self "$species"
done

cat > "${OUT_DIR}/README.txt" <<EOF
All prepared species were analyzed with JCVI self-synteny.
Input root: ${SOURCE_INPUT_ROOT}
Parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}
Species: ${SPECIES[*]}
Raw and screened PNG/PDF dotplots are in ${FIG_DIR}.
EOF

echo "[OK] All-species self-synteny completed."
echo "Summary: ${SUMMARY_DIR}/self_synteny_summary.tsv"
