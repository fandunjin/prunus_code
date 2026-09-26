#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
JCVI_DIR="${KS_ROOT}/jcvi"
KAKS_CONDA_ROOT="${KAKS_CONDA_ROOT:-${CONDA_ROOT}}"
CONDA_SH="${CONDA_SH:-${KAKS_CONDA_ROOT}/etc/profile.d/conda.sh}"
WGD_ENV="${WGD_ENV:-akr}"
KAKS_BIN_DIR="${KAKS_BIN_DIR:-${PROJECT_ROOT}/KaKs_Calculator3.0/src}"
KAKS_CMD="${KAKS_CMD:-${KAKS_BIN_DIR}/KaKs}"
THREADS="${THREADS:-32}"
KAKS_METHOD="${KAKS_METHOD:-YN}"
MIN_AA="${MIN_AA:-50}"
MIN_CODONS="${MIN_CODONS:-50}"
MIN_4D_SITES="${MIN_4D_SITES:-4}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
BOOTSTRAP_SEED="${BOOTSTRAP_SEED:-20260814}"
MIN_CHROM_PAIR_ANCHORS="${MIN_CHROM_PAIR_ANCHORS:-75}"
MIN_CHROM_PAIR_FRACTION="${MIN_CHROM_PAIR_FRACTION:-0.005}"
MIN_REFERENCE_BLOCK_ANCHORS="${MIN_REFERENCE_BLOCK_ANCHORS:-20}"
RUN_LEGACY_BLOCK_REP="${RUN_LEGACY_BLOCK_REP:-0}"

KAKS_SCRIPT="${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py"
BOOTSTRAP_SCRIPT="${ROOT_DIR}/scripts/ks_bootstrap_summary.py"
OUT_ROOT="${KS_ROOT}/kaks"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

mkdir -p "${KS_ROOT}/logs" "${OUT_ROOT}"
[[ -e "${JCVI_DIR}/COMPLETE" ]] || die "Missing completion marker: ${JCVI_DIR}/COMPLETE"
require_file "${JCVI_DIR}/summary/jcvi_syntenic_pairs.tsv"
require_file "${KAKS_SCRIPT}"
require_file "${BOOTSTRAP_SCRIPT}"

set +u
source "${CONDA_SH}"
conda activate "${WGD_ENV}"
set -u
export PATH="${KAKS_BIN_DIR}:${PATH}"

for exe in python mafft pal2nal.pl; do
    command -v "${exe}" >/dev/null 2>&1 || die "${exe} unavailable in ${WGD_ENV}"
done
if [[ ! -x "${KAKS_CMD}" ]]; then
    echo "[ERROR] KaKs_Calculator is not executable: ${KAKS_CMD}" >&2
    echo "[INFO] Contents of ${KAKS_BIN_DIR}:" >&2
    ls -la "${KAKS_BIN_DIR}" >&2 || true
    exit 2
fi

run_analysis() {
    local mode="$1"
    local outdir="${OUT_ROOT}/${mode}_${KAKS_METHOD}"
    rm -f "${outdir}/COMPLETE"
    python "${KAKS_SCRIPT}" \
        --analysis-dir "${KS_ROOT}/stage" \
        --jcvi-dir "${JCVI_DIR}" \
        --outdir "${outdir}" \
        --threads "${THREADS}" \
        --pair-mode "${mode}" \
        --min-chrom-pair-anchors "${MIN_CHROM_PAIR_ANCHORS}" \
        --min-chrom-pair-fraction "${MIN_CHROM_PAIR_FRACTION}" \
        --min-reference-block-anchors "${MIN_REFERENCE_BLOCK_ANCHORS}" \
        --min-aa "${MIN_AA}" \
        --min-codons "${MIN_CODONS}" \
        --min-4d-sites "${MIN_4D_SITES}" \
        --kaks-method "${KAKS_METHOD}" \
        --kaks-cmd "${KAKS_CMD}"

    require_file "${outdir}/summary/all_merged_4dtv_kaks.tsv"
    python "${BOOTSTRAP_SCRIPT}" \
        --input "${outdir}/summary/all_merged_4dtv_kaks.tsv" \
        --out-summary "${outdir}/summary/ks_bootstrap_summary.tsv" \
        --out-contrasts "${outdir}/summary/ks_bootstrap_homeolog_contrasts.tsv" \
        --replicates "${BOOTSTRAPS}" \
        --seed "${BOOTSTRAP_SEED}" \
        --cluster-field block_id \
        --min-ks 0 \
        --max-ks 3 \
        --bin-width 0.02
    touch "${outdir}/COMPLETE"
}

run_analysis collinear_all_pairs

if [[ "${RUN_LEGACY_BLOCK_REP}" == "1" ]]; then
    run_analysis block_representative
fi

{
    echo "conda_root=${KAKS_CONDA_ROOT}"
    echo "conda_sh=${CONDA_SH}"
    echo "kaks_env=${WGD_ENV}"
    echo "kaks_bin_dir=${KAKS_BIN_DIR}"
    echo "kaks_cmd=${KAKS_CMD}"
    echo "python=$(python --version 2>&1)"
    echo "mafft=$(mafft --version 2>&1 | head -n 1)"
    echo "pal2nal=$(command -v pal2nal.pl)"
    echo "KaKs_Calculator=${KAKS_CMD}"
    echo "kaks_method=${KAKS_METHOD}"
    echo "threads=${THREADS}"
    echo "min_aa=${MIN_AA}"
    echo "min_codons=${MIN_CODONS}"
    echo "min_4d_sites=${MIN_4D_SITES}"
    echo "bootstrap_replicates=${BOOTSTRAPS}"
    echo "bootstrap_seed=${BOOTSTRAP_SEED}"
    echo "bootstrap_unit=JCVI_block_id"
    echo "min_chrom_pair_anchors=${MIN_CHROM_PAIR_ANCHORS}"
    echo "min_chrom_pair_fraction=${MIN_CHROM_PAIR_FRACTION}"
    echo "min_reference_block_anchors=${MIN_REFERENCE_BLOCK_ANCHORS}"
} > "${OUT_ROOT}/software_and_parameters.txt"

touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Codon-aware Ks/KaKs analysis completed: ${OUT_ROOT}"
