#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

# Local version of jcvi_pzi_copy1_ordered_self_dotplot.sh. Run from the project
# root on a machine with a working jcvi/diamond Conda environment.
ROOT_DIR="${ROOT_DIR:-$(pwd)}"
HAP_DIR="${HAP_DIR:-${ROOT_DIR}/haploid_rename}"
RENAME_MAP="${RENAME_MAP:-${HAP_DIR}/chromosome_rename_map.tsv}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_ordered_self_jcvi_local}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
JCVI_PYTHON="${JCVI_PYTHON:-}"
THREADS="${THREADS:-8}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"
PREPARE="${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py"
ORDER="5 16 1 2 13 15 11 3 6 12 4 7 8 10 9 14"

die(){ echo "[ERROR] $*" >&2; exit 1; }
req(){ [[ -s "$1" ]] || die "Missing or empty file: $1"; }
count(){ awk 'NF>=2 && $1 !~ /^#/ {n++} END{print n+0}' "$1"; }

CANON="${OUT_DIR}/stage/canonical/pzi_copy1"
INPUT="${OUT_DIR}/inputs"
PAIR="${OUT_DIR}/pairs/pzi_copy1_self"
FIG="${OUT_DIR}/figures"
mkdir -p "$CANON" "$INPUT" "$PAIR" "$FIG" "${OUT_DIR}/summary"
req "$PREPARE"
for h in hapA1 hapB1; do req "${HAP_DIR}/gff/pzi_${h}.gff3"; req "${HAP_DIR}/pep/pzi_${h}.pep.fa"; done
req "$RENAME_MAP"

if [[ "$FORCE" == 1 || ! -s "$CANON/primary.gff3" ]]; then
  { echo '##gff-version 3'; cat "$HAP_DIR/gff/pzi_hapA1.gff3" "$HAP_DIR/gff/pzi_hapB1.gff3"; } > "$CANON/primary.gff3"
fi
if [[ "$FORCE" == 1 || ! -s "$CANON/primary.pep.fa" ]]; then
  cat "$HAP_DIR/pep/pzi_hapA1.pep.fa" "$HAP_DIR/pep/pzi_hapB1.pep.fa" > "$CANON/primary.pep.fa"
fi

# Some Conda activation hooks reference optional variables such as
# LD_LIBRARY_PATH. They are valid during activation but trip `set -u`.
if [[ -n "$JCVI_PYTHON" ]]; then
  [[ -x "$JCVI_PYTHON" ]] || die "JCVI_PYTHON is not executable: $JCVI_PYTHON"
  export PATH="$(dirname "$JCVI_PYTHON"):$PATH"
elif [[ -s "$CONDA_SH" ]]; then
  set +u
  source "$CONDA_SH"
  conda activate "$JCVI_ENV"
  set -u
fi
command -v python >/dev/null || die "python not found; activate jcvi first"
command -v diamond >/dev/null || die "diamond not found in the active environment"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable"
export MPLBACKEND=Agg MPLCONFIGDIR="${OUT_DIR}/mplconfig"; mkdir -p "$MPLCONFIGDIR"

python "$PREPARE" --analysis-dir "${OUT_DIR}/stage" --outdir "$INPUT" \
  --summary "${OUT_DIR}/summary/input_summary.tsv" --labels pzi_copy1
SRC_BED="$INPUT/pzi_copy1/pzi_copy1.bed"; SRC_PEP="$INPUT/pzi_copy1/pzi_copy1.pep"
req "$SRC_BED"; req "$SRC_PEP"

# Rewrite only the BED chromosome order for plotting and JCVI screening.
awk -F '\t' -v OFS='\t' -v order="$ORDER" -v mapfile="$RENAME_MAP" '
 BEGIN{
   n=split(order,a," "); for(i=1;i<=n;i++)r[a[i]]=i
   while ((getline line < mapfile) > 0) {
     if (line ~ /^original_chr32/) continue
     split(line,m,"\t")
     if (m[3] != "" && m[5] != "") rename[m[3]]=m[5]
   }
   close(mapfile)
 }
 NF>=4{
   old=$1; new=rename[old]
   if (new == "") next
   $1=new
   z=new; sub(/^chr/,"",z); z+=0
   if(r[z]) b[r[z],++c[r[z]]]=$0
 }
 END{for(i=1;i<=n;i++)for(j=1;j<=c[i];j++)print b[i,j]}
' "$SRC_BED" > "$PAIR/pzi_copy1.bed.tmp"
mv "$PAIR/pzi_copy1.bed.tmp" "$PAIR/pzi_copy1.bed"
ln -sfn "$SRC_PEP" "$PAIR/pzi_copy1.pep"
cd "$PAIR"
RAW=pzi_copy1.pzi_copy1.anchors; SCREEN=pzi_copy1.pzi_copy1.anchors.new
if [[ "$FORCE" == 1 ]]; then
  rm -f "$RAW" "$SCREEN" *.dotplot.png *.dotplot.pdf \
    "$FIG/pzi_copy1_ordered_self.raw.dotplot.png" \
    "$FIG/pzi_copy1_ordered_self.raw.dotplot.pdf" \
    "$FIG/pzi_copy1_ordered_self.screened.dotplot.png" \
    "$FIG/pzi_copy1_ordered_self.screened.dotplot.pdf"
fi
[[ -s "$RAW" ]] || python -m jcvi.compara.catalog ortholog pzi_copy1 pzi_copy1 --dbtype=prot --align_soft=diamond_blastp --cscore="$CSCORE" --min_size="$MIN_ANCHORS" --cpus="$THREADS" --no_strip_names --no_dotplot
req "$RAW"
[[ -s "$SCREEN" ]] || python -m jcvi.compara.synteny screen --qbed=pzi_copy1.bed --sbed=pzi_copy1.bed --minspan="$MINSPAN" --simple "$RAW" "$SCREEN"
req "$SCREEN"
for x in "$RAW|$FIG/pzi_copy1_ordered_self.raw.dotplot.png|png" "$RAW|$FIG/pzi_copy1_ordered_self.raw.dotplot.pdf|pdf" "$SCREEN|$FIG/pzi_copy1_ordered_self.screened.dotplot.png|png" "$SCREEN|$FIG/pzi_copy1_ordered_self.screened.dotplot.pdf|pdf"; do
 IFS='|' read -r a o f <<< "$x"; [[ -s "$o" ]] || python -m jcvi.graphics.dotplot "$a" --qbed=pzi_copy1.bed --sbed=pzi_copy1.bed --nosort --figsize=12x12 --dpi=300 --format="$f" --nochpf --notex --outfile="$o"; req "$o"
done
printf 'order\traw_anchor_pairs\tscreened_anchor_pairs\n%s\t%s\t%s\n' "$ORDER" "$(count "$RAW")" "$(count "$SCREEN")" > "${OUT_DIR}/summary/summary.tsv"
echo "[OK] Local ordered PZI copy1 dotplots: $FIG"
