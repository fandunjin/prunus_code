#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PAN_LIBRARY="${PAN_LIBRARY:-${PAN_ROOT}/library/frozen/prunus_pan_te_v1.fa}"
LIBRARY_MANIFEST="${LIBRARY_MANIFEST:-${PAN_ROOT}/library/frozen_library_manifest.tsv}"
MASK_ROOT="${MASK_ROOT:-${PAN_ROOT}/mask}"
SELECTION="${SELECTION:-${ROOT_DIR}/0903/TE/family_analysis/subfamily_tree_selection.tsv}"
OUTDIR="${OUTDIR:-${ROOT_DIR}/analysis/prunus_pan_te/followup/subfamily_trees}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_te_subfamily_pipeline.py}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-phylo}"
TOOLS_ENV="${TOOLS_ENV:-soloLTR}"
THREADS="${THREADS:-24}"
MUTATOR_DOMAIN_ROOT="${MUTATOR_DOMAIN_ROOT:-${PAN_ROOT}/mutator_domains}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: missing ${label}: ${path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization"
require_file "${PAN_LIBRARY}" "frozen pan-TE FASTA"
require_file "${LIBRARY_MANIFEST}" "frozen pan-TE manifest"
require_file "${SELECTION}" "subfamily tree selection table"
require_file "${SCRIPT}" "subfamily pipeline script"

source "${CONDA_SH}"
conda activate "${TOOLS_ENV}"
BEDTOOLS_BIN="$(command -v bedtools || true)"
CDHIT_BIN="$(command -v cd-hit-est || true)"
MMSEQS_BIN="$(command -v mmseqs || true)"
TESORTER_BIN="$(command -v TEsorter || true)"
[ -n "${BEDTOOLS_BIN}" ] || { echo "ERROR: bedtools unavailable in ${TOOLS_ENV}" >&2; exit 1; }
[ -n "${TESORTER_BIN}" ] || { echo "ERROR: TEsorter unavailable in ${TOOLS_ENV}" >&2; exit 1; }
if [ -z "${MMSEQS_BIN}" ] && [ -z "${CDHIT_BIN}" ]; then
    echo "ERROR: neither mmseqs nor cd-hit-est is available in ${TOOLS_ENV}" >&2
    exit 1
fi

conda activate "${PYTHON_ENV}"
PYTHON_BIN="$(command -v python)"
MAFFT_BIN="$(command -v mafft || true)"
IQTREE_BIN=""
for candidate in iqtree2 iqtree3 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then
        IQTREE_BIN="$(command -v "${candidate}")"
        break
    fi
done
[ -n "${MAFFT_BIN}" ] || { echo "ERROR: MAFFT unavailable in ${PYTHON_ENV}" >&2; exit 1; }
[ -n "${IQTREE_BIN}" ] || { echo "ERROR: IQ-TREE unavailable in ${PYTHON_ENV}" >&2; exit 1; }

mkdir -p "${OUTDIR}"
args=(
    --selection "${SELECTION}"
    --pan-library "${PAN_LIBRARY}"
    --library-manifest "${LIBRARY_MANIFEST}"
    --mask-root "${MASK_ROOT}"
    --outdir "${OUTDIR}"
    --bedtools "${BEDTOOLS_BIN}"
    --cdhit "${CDHIT_BIN:-cd-hit-est}"
    --tesorter "${TESORTER_BIN}"
    --mafft "${MAFFT_BIN}"
    --iqtree "${IQTREE_BIN}"
    --threads "${THREADS}"
)
if [ -n "${MMSEQS_BIN}" ]; then
    args+=(--mmseqs "${MMSEQS_BIN}")
fi
if [ -d "${MUTATOR_DOMAIN_ROOT}" ]; then
    args+=(--mutator-domain-root "${MUTATOR_DOMAIN_ROOT}")
else
    echo "WARNING: no MUTATOR_DOMAIN_ROOT; DTM families will be reported as missing_mutator_domain." >&2
fi

"${PYTHON_BIN}" "${SCRIPT}" "${args[@]}"
echo "=== TE subfamily tree pipeline completed ==="
echo "Status table: ${OUTDIR}/subfamily_tree_status.tsv"
