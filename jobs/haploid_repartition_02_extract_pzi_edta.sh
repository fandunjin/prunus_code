#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Deprecated PZI-only extractor. Use haploid_repartition_02_extract_edta_all.sh,
# which writes the canonical haploid_repartition/inputs/*/edta tree for all species.

set -euo pipefail

# This job replaces four independent PZI EDTA runs. It subsets the frozen,
# full 2n=32 EDTA annotation after verifying every selected sequence exactly.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/../prunus_old/annotation/repeats/chr32_edta}"
INPUT_DIR="${INPUT_DIR:-${ROOT_DIR}/haploid_repartition/inputs/pzi/genome}"
PZI_MAP="${PZI_MAP:-${ROOT_DIR}/haploid_repartition/maps/pzi_repartition_map.tsv}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/repartition/repeats/pzi}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

SOURCE_FASTA="${SOURCE_FASTA:-${SOURCE_DIR}/pzi_2n32.chromosomes.fa}"
SOURCE_MASKED="${SOURCE_MASKED:-${SOURCE_DIR}/pzi_2n32.chromosomes.fa.mod.MAKER.masked}"
SOURCE_GFF="${SOURCE_GFF:-${SOURCE_DIR}/pzi_2n32.chromosomes.fa.mod.EDTA.TEanno.gff3}"
SOURCE_LIBRARY="${SOURCE_LIBRARY:-${SOURCE_DIR}/pzi_2n32.chromosomes.fa.mod.EDTA.TElib.fa}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/extract_repartitioned_pzi_edta.py}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${SCRIPT}" "EDTA extraction script"
require_file "${SOURCE_FASTA}" "2n=32 EDTA input FASTA"
require_file "${SOURCE_MASKED}" "2n=32 EDTA masked FASTA"
require_file "${SOURCE_GFF}" "2n=32 EDTA TEanno GFF3"
require_file "${SOURCE_LIBRARY}" "2n=32 EDTA library"
require_file "${PZI_MAP}" "current PZI repartition map"
for hap in hapA1 hapA2 hapB1 hapB2; do
    require_file "${INPUT_DIR}/pzi_${hap}.fa" "repartitioned PZI ${hap} FASTA"
done
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable: ${PYTHON_BIN}"

mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"
args=(
    "${SCRIPT}"
    --source-fasta "${SOURCE_FASTA}"
    --source-masked "${SOURCE_MASKED}"
    --source-gff "${SOURCE_GFF}"
    --source-library "${SOURCE_LIBRARY}"
    --pzi-map "${PZI_MAP}"
    --input-genome-dir "${INPUT_DIR}"
    --out-root "${OUT_ROOT}"
)
[[ "${FORCE}" == "1" ]] && args+=(--force)

echo "=== Extract PZI revised-haplotype EDTA annotations ==="
echo "Source EDTA : ${SOURCE_DIR}"
echo "PZI map     : ${PZI_MAP}"
echo "Inputs      : ${INPUT_DIR}"
echo "Outputs     : ${OUT_ROOT}"
echo "Mode        : sequence-verified subset; shared 2n=32 TE library"
"${PYTHON_BIN}" "${args[@]}"

echo "[OK] Do not run haploid_repartition_02_edta.sh for PZI; use this extracted annotation for PZI LAI, Helixer and EviAnn."
