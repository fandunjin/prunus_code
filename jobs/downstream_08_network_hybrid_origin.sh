#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
NETWORK_ENV="${NETWORK_ENV:-network}"
NETWORK_DIR="${ANALYSIS_DIR}/network_hybrid_origin"
TREE_DIR="${TREE_DIR:-${ANALYSIS_DIR}/trees/main}"
GENE_TREES="${GENE_TREES:-${TREE_DIR}/gene_trees/all_gene_trees.tre}"
PHYLONET_NEXUS="${PHYLONET_NEXUS:-${NETWORK_DIR}/phylonet_infernetwork.nex}"
PHYLONET_STATS="${PHYLONET_STATS:-${NETWORK_DIR}/phylonet_gene_tree_filter_stats.tsv}"
PHYLONET_RETICULATIONS="${PHYLONET_RETICULATIONS:-1}"
PHYLONET_OUTGROUPS="${PHYLONET_OUTGROUPS:-fve,rch}"
PHYLONET_FALLBACK_OUTGROUP="${PHYLONET_FALLBACK_OUTGROUP:-fve}"
PHYLONET_MAX_TREES="${PHYLONET_MAX_TREES:-0}"
THREADS="${THREADS:-16}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SPECIES="${SPECIES:-pzi_A2B2,ppe,pmu_refseq,pdu,pcam,pspe,mdom,msyl,pcom,rch,fve}"
PHYLONET_JAR="${PHYLONET_JAR:-${PROJECT_ROOT}/PhyloNet.jar}"
DSUITE_BIN="${DSUITE_BIN:-${PROJECT_ROOT}/Dsuite-master/Build/Dsuite}"

mkdir -p "${NETWORK_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${NETWORK_ENV}"
set -u

if [[ ! -s "${PHYLONET_JAR}" ]]; then
    echo "[ERROR] PHYLONET_JAR is not set to an existing PhyloNet jar: ${PHYLONET_JAR}" >&2
    echo "[ERROR] Set PHYLONET_JAR=/path/to/PhyloNet.jar before submitting downstream_08." >&2
    exit 2
fi

for exe in java run_hyde.py run_hyde_mp.py julia; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in NETWORK_ENV=${NETWORK_ENV}: ${exe}" >&2
        exit 2
    fi
done

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN=python
    else
        echo "[ERROR] Neither ${PYTHON_BIN} nor python is available in NETWORK_ENV=${NETWORK_ENV}" >&2
        exit 2
    fi
fi

python -c "import phyde" >/dev/null 2>&1 || {
    echo "[ERROR] Python package phyde is not importable in NETWORK_ENV=${NETWORK_ENV}" >&2
    exit 2
}

julia -e 'using PhyloNetworks' >/dev/null 2>&1 || {
    echo "[ERROR] Julia package PhyloNetworks is not available in NETWORK_ENV=${NETWORK_ENV}" >&2
    exit 2
}

if [[ ! -x "${DSUITE_BIN}" ]]; then
    echo "[ERROR] DSUITE_BIN is not executable: ${DSUITE_BIN}" >&2
    exit 2
fi

if [ ! -s "${GENE_TREES}" ]; then
    alt_gene_trees="${ANALYSIS_DIR}/TREE/all_gene_trees.tre"
    if [ -s "${alt_gene_trees}" ]; then
        GENE_TREES="${alt_gene_trees}"
    fi
fi

if [ ! -s "${GENE_TREES}" ]; then
    echo "[ERROR] Missing single-copy gene tree collection: ${GENE_TREES}" >&2
    echo "[ERROR] Expected all_gene_trees.tre from downstream_04 under ${TREE_DIR}/gene_trees/." >&2
    exit 2
fi

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_make_phylonet_nexus.py" \
    --gene-trees "${GENE_TREES}" \
    --out-nexus "${PHYLONET_NEXUS}" \
    --stats "${PHYLONET_STATS}" \
    --species "${SPECIES}" \
    --outgroups "${PHYLONET_OUTGROUPS}" \
    --fallback-outgroup "${PHYLONET_FALLBACK_OUTGROUP}" \
    --reticulations "${PHYLONET_RETICULATIONS}" \
    --threads "${THREADS}" \
    --max-trees "${PHYLONET_MAX_TREES}"

if [ ! -s "${PHYLONET_NEXUS}" ]; then
    echo "[ERROR] Failed to generate PhyloNet input: ${PHYLONET_NEXUS}" >&2
    exit 2
fi

java -jar "${PHYLONET_JAR}" "${PHYLONET_NEXUS}" \
    > "${NETWORK_DIR}/phylonet.out" \
    2> "${NETWORK_DIR}/phylonet.err"
