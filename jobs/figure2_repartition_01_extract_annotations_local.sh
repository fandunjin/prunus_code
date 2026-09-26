#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Extract GFF3/CDS/protein files for the revised chromosome packages.
# Run from the project root in WSL. No annotation is re-run.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

require_file() {
    [[ -s "$1" ]] || { echo "[ERROR] Missing or empty input: $1" >&2; exit 2; }
}

require_file "${ROOT_DIR}/haploid_repartition/maps/pzi_repartition_map.tsv"
require_file "${ROOT_DIR}/haploid_repartition/maps/psp_repartition_map.tsv"
require_file "${ROOT_DIR}/haploid_repartition/maps/ppa_repartition_map.tsv"
require_file "${ROOT_DIR}/scripts/extract_repartitioned_haplotype_annotations.py"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || { echo "[ERROR] Python is unavailable: ${PYTHON_BIN}" >&2; exit 2; }

args=("${ROOT_DIR}/scripts/extract_repartitioned_haplotype_annotations.py" --project-root "${ROOT_DIR}")
[[ "${FORCE}" == "1" ]] && args+=(--force)

"${PYTHON_BIN}" "${args[@]}"

echo "[OK] GFF3/CDS/PEP are under haploid_repartition/inputs/{pzi,psp,ppa}/{gff,cds,pep}/"
