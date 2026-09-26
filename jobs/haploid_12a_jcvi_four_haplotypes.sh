#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/circos}"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_pzi_four_haplotype_circos_jcvi.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"

# Relaxed enough to retain fragmented rearrangement blocks (notably C01-C03
# and the provisional C06 group), while retaining a reciprocal-best-hit score
# requirement for protein anchors.
# Apply one common, deliberately permissive filter to all four comparisons
# (A1-A2, B1-B2, A1-B1 and A2-B2).  Keeping these parameters identical makes
# gray and highlighted links directly comparable.
CSCORE="${CSCORE:-0.63}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-12}"
FORCE_JCVI="${FORCE_JCVI:-0}"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"

WORK_DIR="${OUT_DIR}/work"
INPUT_DIR="${WORK_DIR}/jcvi_inputs"
JCVI_DIR="${WORK_DIR}/jcvi"
SUMMARY="${OUT_DIR}/jcvi_block_counts.tsv"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

run_jcvi_pair() {
    local left="$1"
    local right="$2"
    local pair="${left}_${right}"
    local left_label="hap${left}"
    local right_label="hap${right}"
    local pair_dir="${JCVI_DIR}/${pair}"
    local raw="${left_label}.${right_label}.anchors"
    local screened="${left_label}.${right_label}.anchors.screened"
    local search="${left_label}.${right_label}.last"
    local filtered="${search}.filtered"
    local lifted="${left_label}.${right_label}.lifted.anchors"
    local blocks="${pair_dir}/blocks.simple"
    local marker="${pair_dir}/.jcvi.ok"
    local signature="cscore=${CSCORE};min_anchors=${MIN_ANCHORS};minspan=${MINSPAN};package=${PACKAGE_DIGEST}"
    local simple_source=""

    mkdir -p "${pair_dir}"
    ln -sfn "${INPUT_DIR}/${left_label}.bed" "${pair_dir}/${left_label}.bed"
    ln -sfn "${INPUT_DIR}/${left_label}.pep" "${pair_dir}/${left_label}.pep"
    ln -sfn "${INPUT_DIR}/${right_label}.bed" "${pair_dir}/${right_label}.bed"
    ln -sfn "${INPUT_DIR}/${right_label}.pep" "${pair_dir}/${right_label}.pep"

    if [[ "${FORCE_JCVI}" != "1" && -s "${marker}" && -s "${blocks}" ]]; then
        if grep -Fqx -- "${signature}" "${marker}"; then
            echo "[SKIP] Existing JCVI blocks: ${pair}"
            return
        fi
    fi

    cd "${pair_dir}"
    # A changed cscore/min_size invalidates the filtered hits and lifted
    # anchors, but not the expensive all-vs-all DIAMOND search.  Removing only
    # the final .anchors leaves JCVI in an inconsistent cached state, where it
    # skips every stage and never recreates the requested file.
    rm -f "${raw}" "${screened}" "${raw}.simple" "${screened}.simple" \
        "${left_label}.${right_label}.simple" "${filtered}" \
        "${lifted}" "${lifted}.simple" "${blocks}" "${marker}"

    echo
    echo "=== JCVI ${left_label} vs ${right_label} ==="
    python -m jcvi.compara.catalog ortholog \
        "${left_label}" "${right_label}" \
        --dbtype=prot \
        --align_soft=diamond_blastp \
        --cscore="${CSCORE}" \
        --min_size="${MIN_ANCHORS}" \
        --cpus="${THREADS}" \
        --no_strip_names \
        --no_dotplot

    # Some JCVI versions can still accept a stale raw-search cache without
    # rebuilding the final anchors.  Retry once from a completely clean search
    # state only when the expected anchors are still absent.
    if [[ ! -s "${raw}" ]]; then
        echo "[WARN] JCVI cache did not recreate ${raw}; retrying the pair from the sequence search"
        rm -f "${search}" "${filtered}" "${lifted}" "${lifted}.simple" \
            "${raw}" "${raw}.simple"
        python -m jcvi.compara.catalog ortholog \
            "${left_label}" "${right_label}" \
            --dbtype=prot \
            --align_soft=diamond_blastp \
            --cscore="${CSCORE}" \
            --min_size="${MIN_ANCHORS}" \
            --cpus="${THREADS}" \
            --no_strip_names \
            --no_dotplot
    fi
    require_file "${raw}"

    python -m jcvi.compara.synteny screen \
        --qbed="${left_label}.bed" \
        --sbed="${right_label}.bed" \
        --minspan="${MINSPAN}" \
        --simple \
        "${raw}" "${screened}"
    require_file "${screened}"

    for candidate in "${raw}.simple" "${screened}.simple" "${left_label}.${right_label}.simple"; do
        if [[ -s "${candidate}" ]]; then
            simple_source="${candidate}"
            break
        fi
    done
    if [[ -z "${simple_source}" ]]; then
        for candidate in ./*.simple; do
            if [[ -s "${candidate}" ]]; then
                simple_source="${candidate}"
                break
            fi
        done
    fi
    [[ -n "${simple_source}" ]] || die "JCVI did not produce a simple block file for ${pair}"
    ln -sfn "${simple_source}" "${blocks}"
    require_file "${blocks}"
    printf '%s\n' "${signature}" > "${marker}"
    printf '[DONE] %s: ' "${pair}"
    awk 'NF == 6 && $1 !~ /^#/ {n++} END {print n+0, "blocks"}' "${blocks}"
}

mkdir -p "${OUT_DIR}" "${INPUT_DIR}" "${JCVI_DIR}" "${ROOT_DIR}/logs"
require_file "${CONDA_SH}"
require_file "${PREP_SCRIPT}"
require_file "${HAP_RENAME_DIR}/chromosome_rename_map.tsv"
require_file "${PACKAGE_CHECKSUMS}"
for hap in A1 A2 B1 B2; do
    require_file "${HAP_RENAME_DIR}/genome/pzi_hap${hap}.fa"
    require_file "${HAP_RENAME_DIR}/gff/pzi_hap${hap}.gff3"
    require_file "${HAP_RENAME_DIR}/cds/pzi_hap${hap}.cds.fa"
    require_file "${HAP_RENAME_DIR}/pep/pzi_hap${hap}.pep.fa"
done
PACKAGE_DIGEST="$(cksum "${PACKAGE_CHECKSUMS}" | awk '{print $1 ":" $2}')"

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

echo "=== P. zippeliana four-haplotype JCVI preparation ==="
echo "Root       : ${ROOT_DIR}"
echo "Hap package: ${HAP_RENAME_DIR}"
echo "Output     : ${OUT_DIR}"
echo "JCVI env   : ${JCVI_ENV}"
echo "C-score    : ${CSCORE}"
echo "Min anchors: ${MIN_ANCHORS}"
echo "Min span   : ${MINSPAN} genes"
echo "Package ID : ${PACKAGE_DIGEST}"

python "${PREP_SCRIPT}" \
    --root "${ROOT_DIR}" \
    --hap-rename-dir "${HAP_RENAME_DIR}" \
    --out-dir "${OUT_DIR}" \
    --stage-jcvi-inputs

for hap in A1 A2 B1 B2; do
    require_file "${INPUT_DIR}/hap${hap}.bed"
    require_file "${INPUT_DIR}/hap${hap}.pep"
done

# Gray background comparisons and highlighted A-B comparisons.
run_jcvi_pair A1 A2
run_jcvi_pair B1 B2
run_jcvi_pair A1 B1
run_jcvi_pair A2 B2

printf 'pair\tblocks\tfile\n' > "${SUMMARY}"
for pair in A1_A2 B1_B2 A1_B1 A2_B2; do
    blocks="${JCVI_DIR}/${pair}/blocks.simple"
    require_file "${blocks}"
    printf '%s\t' "${pair}" >> "${SUMMARY}"
    awk -v file="${blocks}" 'NF == 6 && $1 !~ /^#/ {n++} END {print n+0 "\t" file}' \
        "${blocks}" >> "${SUMMARY}"
done

echo
echo "=== JCVI preparation complete ==="
echo "Blocks : ${JCVI_DIR}"
echo "Summary: ${SUMMARY}"
