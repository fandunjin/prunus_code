#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Figure 2, stages 2--4: neutral external block discovery and candidate-junction validation.
#
# Design rule: copy1/copy2 are established homolog complements.  This workflow
# never assigns a whole current chromosome to A/B, and an external multi-block
# result is never called a fusion/fission without independent validation.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_block_ancestry}"
THREADS="${THREADS:-12}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CONDA_SH="${CONDA_SH:-}"
JCVI_ENV="${JCVI_ENV:-}"

# Established, neutral PZI inputs.  The same copy1 IDs must be used in every
# external comparison and in the copy1/copy2 and copy1/haplotype validation.
COPY1_INPUT="${COPY1_INPUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_internal_synteny/inputs/copy1}"
COPY2_INPUT="${COPY2_INPUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_internal_synteny/inputs/copy2}"
HAP_INPUT_ROOT="${HAP_INPUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_internal_synteny/inputs}"
INTERNAL_PAIR_ROOT="${INTERNAL_PAIR_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_internal_synteny/pairs}"
COPY_MAP="${COPY_MAP:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv}"
COPY1_FASTA="${COPY1_FASTA:-${ROOT_DIR}/final_genome/pzi_copy1.fa}"

# Existing chromosome-scale Prunus inputs.  PCAM/PPE are filtered again so
# non-chromosomal sequences cannot enter this analysis.  GTR was already
# restricted to Chr01--Chr09 when its neutral input was produced.
GTR_INPUT="${GTR_INPUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison/inputs/gtr}"
PCAM_SOURCE="${PCAM_SOURCE:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs/pcam}"
PPE_SOURCE="${PPE_SOURCE:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs/ppe}"
PCAM_CHROMOSOMES="${PCAM_CHROMOSOMES:-CM114280.1,CM114281.1,CM114282.1,CM114283.1,CM114284.1,CM114285.1,CM114286.1,CM114287.1}"
PPE_CHROMOSOMES="${PPE_CHROMOSOMES:-NC_034009.1,NC_034010.1,NC_034011.1,NC_034012.1,NC_034013.1,NC_034014.1,NC_034015.1,NC_034016.1}"

INPUT_ROOT="${OUT_ROOT}/01_inputs"
PAIR_ROOT="${OUT_ROOT}/02_jcvi_anchors"
BLOCK_ROOT="${OUT_ROOT}/03_syntenic_fragments"
VALIDATION_ROOT="${OUT_ROOT}/04_candidate_validation"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing input: $1"; }

if [[ ! -s "${COPY1_FASTA}.fai" && -s "${COPY1_FASTA}" ]]; then
    command -v samtools >/dev/null 2>&1 || die "samtools is required to index ${COPY1_FASTA}"
    samtools faidx "${COPY1_FASTA}"
fi

if [[ -e "${OUT_ROOT}/COMPLETE" ]]; then
    echo "[OK] Existing neutral block discovery: ${OUT_ROOT}"
    echo "[INFO] To rerun with changed parameters, use a new OUT_ROOT; do not overwrite this completed analysis."
    exit 0
fi

if [[ -n "${CONDA_SH}" ]]; then
    # Optional for WSL.  The HPC wrapper activates the project environment.
    source "${CONDA_SH}"
    [[ -z "${JCVI_ENV}" ]] || conda activate "${JCVI_ENV}"
fi
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable on PATH"
"${PYTHON_BIN}" -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is unavailable in ${PYTHON_BIN}"

for f in \
    "${COPY1_INPUT}/copy1.bed" "${COPY1_INPUT}/copy1.pep" "${COPY1_INPUT}/copy1.cds" \
    "${COPY2_INPUT}/copy2.bed" "${COPY_MAP}" "${COPY1_FASTA}" "${COPY1_FASTA}.fai" \
    "${GTR_INPUT}/gtr.bed" "${GTR_INPUT}/gtr.pep" "${GTR_INPUT}/gtr.cds" \
    "${PCAM_SOURCE}/pcam.bed" "${PCAM_SOURCE}/pcam.pep" "${PCAM_SOURCE}/pcam.cds" \
    "${PPE_SOURCE}/ppe.bed" "${PPE_SOURCE}/ppe.pep" "${PPE_SOURCE}/ppe.cds" \
    "${INTERNAL_PAIR_ROOT}/copy1_vs_copy2/copy1.copy2.anchors.new"; do
    req "${f}"
done
for label in hap01 hap02 hap03 hap04; do
    req "${HAP_INPUT_ROOT}/${label}/${label}.bed"
    req "${INTERNAL_PAIR_ROOT}/copy1_vs_${label}/copy1.${label}.anchors.new"
done

mkdir -p "${INPUT_ROOT}" "${PAIR_ROOT}" "${BLOCK_ROOT}" "${VALIDATION_ROOT}"

# copy1 and GTR are already normalized in the required coordinate/ID systems.
for label in copy1; do
    mkdir -p "${INPUT_ROOT}/${label}"
    for ext in bed pep cds; do
        ln -sfn "${COPY1_INPUT}/${label}.${ext}" "${INPUT_ROOT}/${label}/${label}.${ext}"
    done
done
mkdir -p "${INPUT_ROOT}/gtr"
for ext in bed pep cds; do
    ln -sfn "${GTR_INPUT}/gtr.${ext}" "${INPUT_ROOT}/gtr/gtr.${ext}"
done

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_neutral_external_block_inputs.py" \
    --source pcam "${PCAM_SOURCE}" "${PCAM_CHROMOSOMES}" \
    --source ppe "${PPE_SOURCE}" "${PPE_CHROMOSOMES}" \
    --outdir "${INPUT_ROOT}/external_filtered"
for label in pcam ppe; do
    for ext in bed pep cds; do
        req "${INPUT_ROOT}/external_filtered/${label}/${label}.${ext}"
    done
done

run_pair() {
    local reference="$1"
    local reference_dir="$2"
    local work="${PAIR_ROOT}/copy1_vs_${reference}"
    local raw="copy1.${reference}.anchors"
    local screened="copy1.${reference}.anchors.new"
    mkdir -p "${work}"
    for ext in bed pep cds; do
        ln -sfn "${INPUT_ROOT}/copy1/copy1.${ext}" "${work}/copy1.${ext}"
        ln -sfn "${reference_dir}/${reference}.${ext}" "${work}/${reference}.${ext}"
    done
    pushd "${work}" >/dev/null
    if [[ ! -s "${raw}" ]]; then
        "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog copy1 "${reference}" \
            --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" \
            --min_size="${MIN_ANCHORS}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    fi
    req "${work}/${raw}"
    if [[ ! -s "${screened}" ]]; then
        "${PYTHON_BIN}" -m jcvi.compara.synteny screen \
            --qbed=copy1.bed --sbed="${reference}.bed" --minspan="${MINSPAN}" \
            --simple "${raw}" "${screened}"
    fi
    req "${work}/${screened}"
    popd >/dev/null
}

# GTR is recomputed here rather than reusing an anchor file with different PZI
# IDs.  This is essential for coordinate-consistent downstream validation.
run_pair gtr "${INPUT_ROOT}/gtr"
run_pair pcam "${INPUT_ROOT}/external_filtered/pcam"
run_pair ppe "${INPUT_ROOT}/external_filtered/ppe"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_neutral_syntenic_fragments.py" \
    --query-bed "${INPUT_ROOT}/copy1/copy1.bed" \
    --comparison GTR "${PAIR_ROOT}/copy1_vs_gtr/copy1.gtr.anchors.new" "${INPUT_ROOT}/gtr/gtr.bed" \
    --comparison PCAM "${PAIR_ROOT}/copy1_vs_pcam/copy1.pcam.anchors.new" "${INPUT_ROOT}/external_filtered/pcam/pcam.bed" \
    --comparison PPE "${PAIR_ROOT}/copy1_vs_ppe/copy1.ppe.anchors.new" "${INPUT_ROOT}/external_filtered/ppe/ppe.bed" \
    --outdir "${BLOCK_ROOT}" \
    --min-anchors 10 --max-gap-genes 50 --max-gap-bp 300000 \
    --junction-tolerance-bp 300000 --min-external-support 2

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/validate_neutral_junctions.py" \
    --candidate-junctions "${BLOCK_ROOT}/candidate_junctions.tsv" \
    --query-bed "${INPUT_ROOT}/copy1/copy1.bed" \
    --query-fai "${COPY1_FASTA}.fai" \
    --copy2-bed "${COPY2_INPUT}/copy2.bed" \
    --copy-map "${COPY_MAP}" \
    --copy1-copy2-anchors "${INTERNAL_PAIR_ROOT}/copy1_vs_copy2/copy1.copy2.anchors.new" \
    --hap-evidence hap01 "${INTERNAL_PAIR_ROOT}/copy1_vs_hap01/copy1.hap01.anchors.new" "${HAP_INPUT_ROOT}/hap01/hap01.bed" \
    --hap-evidence hap02 "${INTERNAL_PAIR_ROOT}/copy1_vs_hap02/copy1.hap02.anchors.new" "${HAP_INPUT_ROOT}/hap02/hap02.bed" \
    --hap-evidence hap03 "${INTERNAL_PAIR_ROOT}/copy1_vs_hap03/copy1.hap03.anchors.new" "${HAP_INPUT_ROOT}/hap03/hap03.bed" \
    --hap-evidence hap04 "${INTERNAL_PAIR_ROOT}/copy1_vs_hap04/copy1.hap04.anchors.new" "${HAP_INPUT_ROOT}/hap04/hap04.bed" \
    --outdir "${VALIDATION_ROOT}" \
    --flank-genes 20 --min-side-anchors 11 --min-haplotype-support 2

{
    echo "python=$(${PYTHON_BIN} --version 2>&1)"
    echo "jcvi=$(${PYTHON_BIN} -c 'import jcvi; print(getattr(jcvi, "__version__", "unknown"))' 2>/dev/null || echo unavailable)"
    echo "diamond=$(diamond version 2>&1 | head -n 1)"
    echo "cscore=${CSCORE}"
    echo "min_size=${MIN_ANCHORS}"
    echo "minspan=${MINSPAN}"
    echo "external_references=GTR(Chr01-Chr09),PCAM(8 chromosomes),PPE(8 chromosomes)"
    echo "fragment_discovery=min_anchors:10,max_gap_genes:50,max_gap_bp:300000,junction_tolerance_bp:300000,min_external_support:2"
    echo "validation=20 genes each side; >10 conserved anchors on each side; 100-kb review windows"
    echo "claim_boundary=external multiblock results are candidates only; no A/B, LF/MF, fusion/fission, or parental-origin call"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Neutral external block discovery completed: ${OUT_ROOT}"
echo "[NEXT] Review ${VALIDATION_ROOT}/junction_validation_summary.tsv and validate retained candidates with correctly referenced HiFi/Hi-C evidence."
