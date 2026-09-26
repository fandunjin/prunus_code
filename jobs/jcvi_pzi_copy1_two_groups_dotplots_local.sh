#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

# Local JCVI dotplots for two selected PZI copy1 chromosome groups.
# Group 1: 5,1,2,3,4,8,9
# Group 2: 16,13,15,11,6,12,7,10,14
# Produces group1 self, group2 self and group1-vs-group2 dotplots.

ROOT_DIR="${ROOT_DIR:-$(pwd)}"
HAP_DIR="${HAP_DIR:-${ROOT_DIR}/haploid_rename}"
RENAME_MAP="${RENAME_MAP:-${HAP_DIR}/chromosome_rename_map.tsv}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_copy1_two_groups_jcvi_local}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
JCVI_PYTHON="${JCVI_PYTHON:-}"
THREADS="${THREADS:-8}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"
PREPARE="${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py"
GROUP1="${GROUP1:-1 2 5 6 7 10 12 14}"
GROUP2="${GROUP2:-13 15 11 16 3 4 8 9}"

die(){ echo "[ERROR] $*" >&2; exit 1; }
req(){ [[ -s "$1" ]] || die "Missing or empty file: $1"; }
count(){ awk 'NF>=2 && $1 !~ /^#/ {n++} END{print n+0}' "$1"; }

CANON="${OUT_DIR}/stage/canonical/pzi_copy1"
INPUT="${OUT_DIR}/inputs"
PAIR_ROOT="${OUT_DIR}/pairs"
FIG="${OUT_DIR}/figures"
SUMMARY="${OUT_DIR}/summary"
mkdir -p "$CANON" "$INPUT" "$PAIR_ROOT" "$FIG" "$SUMMARY"
req "$PREPARE"; req "$RENAME_MAP"
for h in hapA1 hapB1; do
  req "${HAP_DIR}/gff/pzi_${h}.gff3"
  req "${HAP_DIR}/pep/pzi_${h}.pep.fa"
done

if [[ "$FORCE" == 1 || ! -s "$CANON/primary.gff3" ]]; then
  { printf '##gff-version 3\n'; cat "$HAP_DIR/gff/pzi_hapA1.gff3" "$HAP_DIR/gff/pzi_hapB1.gff3"; } > "$CANON/primary.gff3"
fi
if [[ "$FORCE" == 1 || ! -s "$CANON/primary.pep.fa" ]]; then
  cat "$HAP_DIR/pep/pzi_hapA1.pep.fa" "$HAP_DIR/pep/pzi_hapB1.pep.fa" > "$CANON/primary.pep.fa"
fi

if [[ -n "$JCVI_PYTHON" ]]; then
  [[ -x "$JCVI_PYTHON" ]] || die "JCVI_PYTHON is not executable: $JCVI_PYTHON"
  export PATH="$(dirname "$JCVI_PYTHON"):$PATH"
elif [[ -s "$CONDA_SH" ]]; then
  set +u; source "$CONDA_SH"; conda activate "$JCVI_ENV"; set -u
fi
command -v python >/dev/null || die "python not found"
command -v diamond >/dev/null || die "diamond not found"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable"
export MPLBACKEND=Agg MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"

python "$PREPARE" --analysis-dir "${OUT_DIR}/stage" --outdir "$INPUT" \
  --summary "$SUMMARY/input_summary.tsv" --labels pzi_copy1
SRC_BED="$INPUT/pzi_copy1/pzi_copy1.bed"
SRC_PEP="$INPUT/pzi_copy1/pzi_copy1.pep"
req "$SRC_BED"; req "$SRC_PEP"

# Convert haplotype chromosome IDs to copy1 chr01--chr16 IDs and retain all
# genes. This is the common source BED for the two extracted groups.
awk -F '\t' -v OFS='\t' -v mapfile="$RENAME_MAP" '
 BEGIN {
   while ((getline line < mapfile) > 0) {
     if (line ~ /^original_chr32/) continue
     split(line,m,"\t")
     if (m[3] != "" && m[5] != "") rename[m[3]]=m[5]
   }
   close(mapfile)
 }
 NF >= 4 {
   new=rename[$1]
   if (new == "") next
   $1=new
   print
 }
' "$SRC_BED" > "$OUT_DIR/copy1_all.bed"
req "$OUT_DIR/copy1_all.bed"

extract_group() {
  local name="$1" order="$2"
  local bed="$INPUT/${name}/${name}.bed"
  local pep="$INPUT/${name}/${name}.pep"
  mkdir -p "$INPUT/$name"
  awk -F '\t' -v OFS='\t' -v order="$order" '
    BEGIN{n=split(order,a," "); for(i=1;i<=n;i++)rank["chr" sprintf("%02d",a[i])]=i}
    NF>=4 && rank[$1] {b[rank[$1],++c[rank[$1]]]=$0}
    END{for(i=1;i<=n;i++)for(j=1;j<=c[i];j++)print b[i,j]}
  ' "$OUT_DIR/copy1_all.bed" > "$bed"
  req "$bed"
  awk -F '\t' 'NF>=4 {print $4}' "$bed" > "$INPUT/${name}/ids.txt"
  awk -v ids="$INPUT/${name}/ids.txt" '
    BEGIN{while((getline x < ids)>0)keep[x]=1; close(ids)}
    /^>/ {
      if(out!="" && keep[id]) printf "%s", out
      id=$1; sub(/^>/,"",id); split(id,a,/[^[:alnum:]_.:-]/); id=a[1]
      out=(keep[id] ? $0 "\n" : ""); next
    }
    {if(out!="") out=out $0 "\n"}
    END{if(out!="" && keep[id]) printf "%s", out}
  ' "$SRC_PEP" > "$pep"
  req "$pep"
}

extract_group pzi_copy1_group1 "$GROUP1"
extract_group pzi_copy1_group2 "$GROUP2"

run_pair() {
  local query="$1" ref="$2" tag="$3"
  local pair="${PAIR_ROOT}/${tag}"
  local qbed="${query}.bed" rbed="${ref}.bed"
  local raw="${query}.${ref}.anchors" screen="${query}.${ref}.anchors.new"
  local rawpng="$FIG/${tag}.raw.dotplot.png" rawpdf="$FIG/${tag}.raw.dotplot.pdf"
  local scrpng="$FIG/${tag}.screened.dotplot.png" scrpdf="$FIG/${tag}.screened.dotplot.pdf"
  mkdir -p "$pair"
  ln -sfn "$INPUT/$query/$qbed" "$pair/$qbed"
  ln -sfn "$INPUT/$ref/$rbed" "$pair/$rbed"
  ln -sfn "$INPUT/$query/${query}.pep" "$pair/${query}.pep"
  ln -sfn "$INPUT/$ref/${ref}.pep" "$pair/${ref}.pep"
  if [[ "$FORCE" == 1 ]]; then rm -f "$pair/$raw" "$pair/$screen" "$pair"/*.simple "$rawpng" "$rawpdf" "$scrpng" "$scrpdf"; fi
  (
    cd "$pair"
    [[ -s "$raw" ]] || python -m jcvi.compara.catalog ortholog "$query" "$ref" \
      --dbtype=prot --align_soft=diamond_blastp --cscore="$CSCORE" \
      --min_size="$MIN_ANCHORS" --cpus="$THREADS" --no_strip_names --no_dotplot
    req "$raw"
    [[ -s "$screen" ]] || python -m jcvi.compara.synteny screen \
      --qbed="$qbed" --sbed="$rbed" --minspan="$MINSPAN" --simple "$raw" "$screen"
    req "$screen"
    for spec in "$raw|$rawpng|png" "$raw|$rawpdf|pdf" "$screen|$scrpng|png" "$screen|$scrpdf|pdf"; do
      IFS='|' read -r a o f <<< "$spec"
      [[ -s "$o" ]] || python -m jcvi.graphics.dotplot "$a" --qbed="$qbed" --sbed="$rbed" \
        --nosort --figsize=12x12 --dpi=300 --format="$f" --nochpf --notex --outfile="$o"
      req "$o"
    done
  )
  printf '%s\t%s\t%s\t%s\t%s\n' "$tag" "$(count "$pair/$raw")" "$(count "$pair/$screen")" "$rawpng" "$scrpng" >> "$SUMMARY/pairwise_summary.tsv"
}

printf 'comparison\traw_anchor_pairs\tscreened_anchor_pairs\traw_dotplot\tscreened_dotplot\n' > "$SUMMARY/pairwise_summary.tsv"
run_pair pzi_copy1_group1 pzi_copy1_group1 group1_self
run_pair pzi_copy1_group2 pzi_copy1_group2 group2_self
run_pair pzi_copy1_group1 pzi_copy1_group2 group1_vs_group2

cat > "$OUT_DIR/README.txt" <<EOF
PZI copy1 selected chromosome-group JCVI dotplots.
Group 1 order: ${GROUP1}
Group 2 order: ${GROUP2}
Comparisons: group1 self, group2 self, group1 versus group2.
Parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}.
EOF
echo "[OK] Three group dotplots written to $FIG"
