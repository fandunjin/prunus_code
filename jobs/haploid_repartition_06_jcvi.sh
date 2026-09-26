#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/circos/repartitioned}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_jcvi_inputs.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.63}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-12}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
req "${SCRIPT}" "JCVI input preparation script"
req "${CONDA_SH}" "conda initialization"
mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python "${SCRIPT}" --root "${ROOT_DIR}" --out-root "${OUT_ROOT}"

INPUT_DIR="${OUT_ROOT}/work/jcvi_inputs"
JCVI_DIR="${OUT_ROOT}/work/jcvi"
mkdir -p "${JCVI_DIR}"

run_pair() {
    local left="$1" right="$2" pair="${left}_${right}"
    local dir="${JCVI_DIR}/${pair}"
    local raw="${left}.${right}.anchors"
    local screened="${raw}.screened"
    local blocks="${dir}/blocks.simple"
    mkdir -p "${dir}"
    ln -sfn "${INPUT_DIR}/${left}.bed" "${dir}/${left}.bed"
    ln -sfn "${INPUT_DIR}/${left}.pep" "${dir}/${left}.pep"
    ln -sfn "${INPUT_DIR}/${right}.bed" "${dir}/${right}.bed"
    ln -sfn "${INPUT_DIR}/${right}.pep" "${dir}/${right}.pep"
    if [[ "${FORCE}" == "1" ]]; then
        rm -f "${dir}/${raw}" "${dir}/${screened}" "${blocks}" "${dir}"/*.last "${dir}"/*.filtered "${dir}"/*.lifted.anchors
    fi
    cd "${dir}"
    if [[ ! -s "${raw}" ]]; then
        python -m jcvi.compara.catalog ortholog "${left}" "${right}" --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" --min_size="${MIN_ANCHORS}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    fi
    req "${dir}/${raw}" "JCVI raw anchors ${pair}"
    python -m jcvi.compara.synteny screen --qbed="${left}.bed" --sbed="${right}.bed" --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
    req "${dir}/${screened}" "JCVI screened anchors ${pair}"
    local source=""
    for candidate in "${raw}.simple" "${screened}.simple" "${left}.${right}.simple"; do
        if [[ -s "${candidate}" ]]; then source="${candidate}"; break; fi
    done
    [[ -n "${source}" ]] || die "No JCVI simple block output for ${pair}"
    ln -sfn "${source}" "${blocks}"
    echo "[DONE] ${pair}: $(awk 'NF==6 && $1 !~ /^#/ {n++} END {print n+0}' "${blocks}") blocks"
}

# All filters are identical across within-haplotype and cross-species pairs.
for pair in \
    "pzi_A1 pzi_A2" "pzi_B1 pzi_B2" "pzi_A1 pzi_B1" "pzi_A2 pzi_B2" \
    "ppa_A ppa_B" \
    "pzi_A1 psp_A1" "pzi_B1 psp_B1" "pzi_A1 ppa_A" "pzi_B1 ppa_B" \
    "psp_A1 ppa_A" "psp_B1 ppa_B"; do
    read -r left right <<< "${pair}"
    run_pair "${left}" "${right}"
done

printf 'pair\tblocks\tfile\n' > "${OUT_ROOT}/jcvi_block_summary.tsv"
for dir in "${JCVI_DIR}"/*; do
    [[ -s "${dir}/blocks.simple" ]] || continue
    pair="$(basename "${dir}")"
    printf '%s\t%s\t%s\n' "${pair}" "$(awk 'NF==6 && $1 !~ /^#/ {n++} END {print n+0}' "${dir}/blocks.simple")" "${dir}/blocks.simple" >> "${OUT_ROOT}/jcvi_block_summary.tsv"
done
echo "[OK] JCVI blocks: ${JCVI_DIR}"
