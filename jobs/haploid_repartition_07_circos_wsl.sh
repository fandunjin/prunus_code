#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL replacement for the native Circos rendering job.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/circos/repartitioned}"
JCVI_ROOT="${JCVI_ROOT:-${OUT_ROOT}/work/jcvi}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_native_circos.py"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
req "${SCRIPT}" "native Circos preparation script"
req "${OUT_ROOT}/work/jcvi_inputs/pzi_A1.bed" "revised PZI JCVI BED"
command -v python >/dev/null 2>&1 || die "python is unavailable in the active WSL environment"
command -v circos >/dev/null 2>&1 || die "circos is unavailable; run haploid_repartition_00_wsl_setup.sh"

if [[ "${FORCE}" == "1" || ! -s "${OUT_ROOT}/pzi/conf/circos.conf" ]]; then
    python "${SCRIPT}" --root "${ROOT_DIR}" --out-root "${OUT_ROOT}" --jcvi-root "${JCVI_ROOT}"
else
    echo "[SKIP] Existing native Circos inputs: ${OUT_ROOT}"
fi

for species in pzi psp ppa; do
    req "${OUT_ROOT}/${species}/conf/circos.conf" "${species} Circos configuration"
    req "${OUT_ROOT}/${species}/data/karyotype.txt" "${species} Circos karyotype"
    req "${OUT_ROOT}/${species}/data/links_highlight.txt" "${species} Circos links"
    cd "${OUT_ROOT}/${species}"
    rm -f "figures/${species}_repartition_circos.png" "figures/${species}_repartition_circos.svg"
    circos -conf conf/circos.conf
    req "figures/${species}_repartition_circos.png" "${species} Circos PNG"
    req "figures/${species}_repartition_circos.svg" "${species} Circos SVG"
    grep -q '</svg>' "figures/${species}_repartition_circos.svg" || die "${species} SVG is incomplete"
    echo "[OK] ${species}: ${OUT_ROOT}/${species}/figures/${species}_repartition_circos.png"
done
