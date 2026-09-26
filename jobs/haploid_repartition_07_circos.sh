#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/circos/repartitioned}"
JCVI_ROOT="${JCVI_ROOT:-${OUT_ROOT}/work/jcvi}"
JCVI_INPUT_ROOT="${JCVI_INPUT_ROOT:-${OUT_ROOT}/work/jcvi_inputs}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_native_circos.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
CIRCOS_ENV="${CIRCOS_ENV:-circos}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
req "${SCRIPT}" "Circos preparation script"
req "${CONDA_SH}" "conda initialization"
req "${JCVI_INPUT_ROOT}/pzi_A1.bed" "staged PZI A1 BED"
req "${ROOT_DIR}/haploid_repartition/inputs/pzi/genome/pzi_hapA1.fa" "revised PZI A1 FASTA"

mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${CIRCOS_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${CIRCOS_ENV}"
command -v circos >/dev/null 2>&1 || die "circos is unavailable in ${CIRCOS_ENV}"

if [[ "${FORCE}" == "1" || ! -s "${OUT_ROOT}/pzi/conf/circos.conf" ]]; then
    python "${SCRIPT}" \
        --root "${ROOT_DIR}" \
        --out-root "${OUT_ROOT}" \
        --jcvi-root "${JCVI_ROOT}"
else
    echo "[SKIP] Existing Circos inputs: ${OUT_ROOT}"
fi

for species in pzi psp ppa; do
    req "${OUT_ROOT}/${species}/conf/circos.conf" "${species} Circos configuration"
    req "${OUT_ROOT}/${species}/data/karyotype.txt" "${species} Circos karyotype"
    req "${OUT_ROOT}/${species}/data/links_highlight.txt" "${species} Circos highlight links"
    cd "${OUT_ROOT}/${species}"
    rm -f figures/${species}_repartition_circos.png figures/${species}_repartition_circos.svg
    circos -conf conf/circos.conf
    req "${OUT_ROOT}/${species}/figures/${species}_repartition_circos.png" "${species} Circos PNG"
    req "${OUT_ROOT}/${species}/figures/${species}_repartition_circos.svg" "${species} Circos SVG"
    grep -q '</svg>' "${OUT_ROOT}/${species}/figures/${species}_repartition_circos.svg" || die "${species} Circos SVG is incomplete"
    echo "[OK] Circos figure: ${OUT_ROOT}/${species}/figures/${species}_repartition_circos.png"
done
