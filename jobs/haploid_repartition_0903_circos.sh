#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT_ROOT}/0903}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/haploid_repartition/circos}"
JCVI_ROOT="${JCVI_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work}"
BED_ROOT="${BED_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks/inputs}"
TRACK_ROOT="${TRACK_ROOT:-${PROJECT_ROOT}/downstream_analysis/circos/repartitioned}"
SCRIPT="${SCRIPT:-${PROJECT_ROOT}/scripts/prepare_0903_ab_circos.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
CIRCOS_ENV="${CIRCOS_ENV:-circos}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

req "${SCRIPT}" "CircOS preparation script"
req "${CONDA_SH}" "conda initialization"
for input in \
    "${ROOT_DIR}/haploid_repartition/inputs/pzi/genome/pzi_hapA1.fa" \
    "${ROOT_DIR}/haploid_repartition/inputs/psp/genome/psp_hapA1.fa" \
    "${ROOT_DIR}/haploid_repartition/inputs/ppa/genome/ppa_hapA.fa" \
    "${JCVI_ROOT}/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors.simple" \
    "${JCVI_ROOT}/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors.simple" \
    "${JCVI_ROOT}/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors.simple"; do
    req "${input}" "required CircOS input"
done

mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${CIRCOS_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${CIRCOS_ENV}"
command -v circos >/dev/null 2>&1 || die "circos is unavailable in ${CIRCOS_ENV}"

CIRCOS_BIN="$(command -v circos)"
CIRCOS_ETC="${CIRCOS_ETC:-}"
if [[ -z "${CIRCOS_ETC}" ]]; then
    for candidate in \
        "$(dirname "${CIRCOS_BIN}")/../etc" \
        "$(dirname "${CIRCOS_BIN}")/etc" \
        "/etc/circos"; do
        if [[ -s "${candidate}/colors_fonts_patterns.conf" ]]; then
            CIRCOS_ETC="$(cd "${candidate}" && pwd)"
            break
        fi
    done
fi
[[ -n "${CIRCOS_ETC}" ]] || die "Cannot locate CircOS etc directory"

if [[ "${FORCE}" == "1" || ! -s "${OUT_ROOT}/figure_manifest.tsv" ]]; then
    python "${SCRIPT}" \
        --root "${ROOT_DIR}" \
        --out-root "${OUT_ROOT}" \
        --jcvi-root "${JCVI_ROOT}" \
        --bed-root "${BED_ROOT}" \
        --track-root "${TRACK_ROOT}"
else
    echo "[SKIP] Existing prepared inputs: ${OUT_ROOT} (set FORCE=1 to refresh)"
fi

for species in pzi psp ppa; do
    species_root="${OUT_ROOT}/${species}"
    req "${species_root}/conf/circos.conf" "${species} CircOS configuration"
    req "${species_root}/data/karyotype.txt" "${species} CircOS karyotype"
    req "${species_root}/data/links_ab.txt" "${species} A-B links"
    mkdir -p "${species_root}/conf"
    ln -sfn "${CIRCOS_ETC}" "${species_root}/conf/etc"
    cd "${species_root}"
    rm -f "figures/${species}_ab_circos.png" "figures/${species}_ab_circos.svg"
    circos -conf conf/circos.conf
    req "${species_root}/figures/${species}_ab_circos.png" "${species} CircOS PNG"
    req "${species_root}/figures/${species}_ab_circos.svg" "${species} CircOS SVG"
    grep -q '</svg>' "${species_root}/figures/${species}_ab_circos.svg" || die "${species} SVG is incomplete"
    req "${species_root}/figures/${species}_ab_circos.pdf" "${species} CircOS PDF"
    echo "[OK] ${species}: ${species_root}/figures/${species}_ab_circos.png"
done

echo "[OK] PZI/PSP/PPA A-B CircOS figures: ${OUT_ROOT}"
