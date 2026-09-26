#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Neutral copy1-vs-GTR synteny and K=8/K=9-block/K=9+7 comparison.
# This job never uses PZI A/B labels or the old subgenome assignment table.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison}"
THREADS="${THREADS:-12}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CONDA_SH="${CONDA_SH:-}"
JCVI_ENV="${JCVI_ENV:-}"
FORCE="${FORCE:-0}"

PZI_GFF="${PZI_GFF:-${ROOT_DIR}/annotation/genes/copy1_helixer_edta_masked/postprocess/pzi_copy1.edta_masked.helixer.genes.gff3}"
PZI_PEP="${PZI_PEP:-${ROOT_DIR}/annotation/genes/copy1_helixer_edta_masked/postprocess/pzi_copy1.edta_masked.helixer.protein.fa}"
PZI_CDS="${PZI_CDS:-${ROOT_DIR}/annotation/genes/copy1_helixer_edta_masked/postprocess/pzi_copy1.edta_masked.helixer.cds.fa}"
PZI_GENOME="${PZI_GENOME:-${ROOT_DIR}/final_genome/pzi_copy1.fa}"
GTR_INPUT="${GTR_INPUT:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs/gtr}"

STAGE="${OUT_ROOT}/stage"
INPUTS="${OUT_ROOT}/inputs"
PAIR="${OUT_ROOT}/pair/pzi_copy1_vs_gtr"
RESULTS="${OUT_ROOT}/model_results"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing input: $1"; }

for f in "${PZI_GFF}" "${PZI_PEP}" "${PZI_CDS}" "${PZI_GENOME}" \
         "${GTR_INPUT}/gtr.bed" "${GTR_INPUT}/gtr.pep" "${GTR_INPUT}/gtr.cds"; do
    req "${f}"
done
if [[ -e "${OUT_ROOT}/COMPLETE" && "${FORCE}" != "1" ]]; then
    echo "[OK] Existing neutral copy1/GTR model comparison: ${OUT_ROOT}"
    exit 0
fi

mkdir -p "${STAGE}/canonical/pzi_copy1" "${INPUTS}" "${PAIR}" "${RESULTS}"
ln -sfn "${PZI_GENOME}" "${STAGE}/canonical/pzi_copy1/genome.fa"
ln -sfn "${PZI_GFF}" "${STAGE}/canonical/pzi_copy1/primary.gff3"
ln -sfn "${PZI_PEP}" "${STAGE}/canonical/pzi_copy1/primary.pep.fa"
ln -sfn "${PZI_CDS}" "${STAGE}/canonical/pzi_copy1/primary.cds.fa"

if [[ -n "${CONDA_SH}" ]]; then
    source "${CONDA_SH}"
    [[ -z "${JCVI_ENV}" ]] || conda activate "${JCVI_ENV}"
fi

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py" \
    --analysis-dir "${STAGE}" \
    --outdir "${INPUTS}" \
    --summary "${OUT_ROOT}/input_summary.tsv" \
    --labels pzi_copy1

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_neutral_copy1_gtr_inputs.py" \
    --gtr-input "${GTR_INPUT}" \
    --outdir "${INPUTS}/gtr"

for f in "${INPUTS}/pzi_copy1/pzi_copy1.bed" "${INPUTS}/pzi_copy1/pzi_copy1.pep" \
         "${INPUTS}/gtr/gtr.bed" "${INPUTS}/gtr/gtr.pep"; do req "${f}"; done

for label in pzi_copy1 gtr; do
    for ext in bed pep cds; do
        src="${INPUTS}/${label}/${label}.${ext}"
        [[ -s "${src}" ]] || die "Normalized ${label}.${ext} is empty"
        ln -sfn "${src}" "${PAIR}/${label}.${ext}"
    done
done

if [[ "${FORCE}" == "1" ]]; then
    rm -f "${PAIR}/pzi_copy1.gtr.anchors" "${PAIR}/pzi_copy1.gtr.anchors.new" \
          "${PAIR}/pzi_copy1.gtr.anchors.simple" "${PAIR}/pzi_copy1.gtr.anchors.new.simple"
fi
pushd "${PAIR}" >/dev/null
if [[ ! -s pzi_copy1.gtr.anchors ]]; then
    "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog pzi_copy1 gtr \
        --dbtype=prot --align_soft=diamond_blastp \
        --cscore="${CSCORE}" --min_size="${MIN_ANCHORS}" \
        --cpus="${THREADS}" --no_strip_names --no_dotplot
fi
req "${PAIR}/pzi_copy1.gtr.anchors"
if [[ ! -s pzi_copy1.gtr.anchors.new ]]; then
    "${PYTHON_BIN}" -m jcvi.compara.synteny screen \
        --qbed=pzi_copy1.bed --sbed=gtr.bed --minspan="${MINSPAN}" \
        --simple pzi_copy1.gtr.anchors pzi_copy1.gtr.anchors.new
fi
req "${PAIR}/pzi_copy1.gtr.anchors.new"
popd >/dev/null

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/compare_neutral_copy1_gtr_models.py" \
    --anchors "${PAIR}/pzi_copy1.gtr.anchors.new" \
    --qbed "${INPUTS}/pzi_copy1/pzi_copy1.bed" \
    --gtr-bed "${INPUTS}/gtr/gtr.bed" \
    --outdir "${RESULTS}"

{
    echo "root=${ROOT_DIR}"
    echo "python=$(${PYTHON_BIN} --version 2>&1)"
    echo "jcvi=$(${PYTHON_BIN} -c 'import jcvi; print(getattr(jcvi, "__version__", "unknown"))' 2>/dev/null || echo unavailable)"
    echo "cscore=${CSCORE}"
    echo "min_size=${MIN_ANCHORS}"
    echo "minspan=${MINSPAN}"
    echo "gtr_chromosomes=Chr01-Chr09 only"
    echo "pzi_labels=neutral copy1 chr01-chr16; no A/B labels"
    echo "comparison_score=100*(0.6*anchor_retention+0.2*order_continuity+0.2*structural_coherence)-2*event_burden-0.5*extra_assignments"
    echo "interpretation=heuristic comparative score; not a formal likelihood"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Neutral copy1/GTR model comparison completed: ${OUT_ROOT}"
echo "[OK] Model table: ${RESULTS}/model_comparison.tsv"
