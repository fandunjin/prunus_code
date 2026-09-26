#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/figure2/01_ppa_pzi_copy1_minimap}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/data/ppa/haploid_rename}"
PPA_GENOME="${PPA_GENOME:-${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
PPA_GFF="${PPA_GFF:-${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.genes.gff3}"
PPA_CDS="${PPA_CDS:-${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.cds.fa}"
PPA_PEP="${PPA_PEP:-${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.protein.fa}"
PRIMARY_MAP="${PRIMARY_MAP:-${MINIMAP_DIR}/ppa_primary_chr16.rename_map.tsv}"
ASSIGNMENT="${ASSIGNMENT:-${MINIMAP_DIR}/ppa_to_pzi_copy1_assignment.tsv}"
BEST_HITS="${BEST_HITS:-${MINIMAP_DIR}/ppa_vs_pzi_copy1_chr16.asm5.best_pzi_copy1_hit_per_ppa_superscaffold.tsv}"
PZI_NAMING_MAP="${PZI_NAMING_MAP:-${ROOT_DIR}/final_haplotype_genomes/pzi_chr32_copy_haplotype_naming.tsv}"
ASSIGNMENT_SCRIPT="${ASSIGNMENT_SCRIPT:-${ROOT_DIR}/scripts/build_ppa_assignment_from_pzi_copy1.py}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/partition_ppa_subgenomes.py}"
FORCE="${FORCE:-0}"
ASSIGNMENT_STATUS="${ASSIGNMENT_STATUS:-approved}"

for path in "${PPA_GENOME}" "${PPA_GFF}" "${PPA_CDS}" "${PPA_PEP}" "${PRIMARY_MAP}" "${BEST_HITS}" "${PZI_NAMING_MAP}" "${ASSIGNMENT_SCRIPT}" "${SCRIPT}"; do
    [[ -s "${path}" ]] || { echo "[ERROR] Missing input: ${path}" >&2; exit 2; }
done
command -v python3 >/dev/null 2>&1 && PYTHON_BIN=python3 || PYTHON_BIN=python

[[ "${ASSIGNMENT_STATUS}" == "approved" || "${ASSIGNMENT_STATUS}" == "pending" ]] || { echo "[ERROR] ASSIGNMENT_STATUS must be approved or pending" >&2; exit 2; }
"${PYTHON_BIN}" "${ASSIGNMENT_SCRIPT}" \
    --best-hits "${BEST_HITS}" \
    --pzi-naming-map "${PZI_NAMING_MAP}" \
    --out "${ASSIGNMENT}" \
    --review-status "${ASSIGNMENT_STATUS}"

args=(--genome "${PPA_GENOME}" --gff "${PPA_GFF}" --cds "${PPA_CDS}" --pep "${PPA_PEP}" --primary-map "${PRIMARY_MAP}" --assignment "${ASSIGNMENT}" --out-dir "${OUT_DIR}")
[[ "${FORCE}" == "1" ]] && args+=(--force)
"${PYTHON_BIN}" "${SCRIPT}" "${args[@]}"

echo "[OK] PPA A/B package: ${OUT_DIR}"
