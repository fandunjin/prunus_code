#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PAML_ENV="${PAML_ENV:-paml}"
IQTREE_ENV="${IQTREE_ENV:-orthofinder}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-16}"
TREE_DIR="${ANALYSIS_DIR}/trees/main"
DATE_DIR="${ANALYSIS_DIR}/dating/mcmctree_4d"
IQ2MC_DIR="${DATE_DIR}/iqtree_mcmctree"
IQ2MC_PREFIX_NAME="${IQ2MC_PREFIX_NAME:-iq2mc_4d}"
IQ2MC_PREFIX="${IQ2MC_DIR}/${IQ2MC_PREFIX_NAME}"
CHAIN_ID="${CHAIN_ID:-chain1}"
CHAIN_DIR="${DATE_DIR}/${CHAIN_ID}"
MCMC_BURNIN="${MCMC_BURNIN:-200000}"
MCMC_SAMPFREQ="${MCMC_SAMPFREQ:-50}"
MCMC_NSAMPLE="${MCMC_NSAMPLE:-50000}"
MCMC_BDS="${MCMC_BDS:-1,1,0.5}"
MCMC_CLOCK="${MCMC_CLOCK:-IND}"
IQTREE_DATING_MODEL="${IQTREE_DATING_MODEL:-MFP}"

mkdir -p "${DATE_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${PAML_ENV}"
set -u

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN=python
    else
        echo "[ERROR] Neither ${PYTHON_BIN} nor python is available in PAML_ENV=${PAML_ENV}" >&2
        exit 2
    fi
fi

for exe in seqmagick; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in PAML_ENV=${PAML_ENV}: ${exe}" >&2
        exit 2
    fi
done

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_extract_4dtv.py" \
    --alignment-dir "${TREE_DIR}/codon_alignments" \
    --pattern "*.fa" \
    --out-prefix "${DATE_DIR}/single_copy"

seqmagick convert "${DATE_DIR}/single_copy.4d_sites.fa" "${DATE_DIR}/single_copy.4d_sites.phy"

cp -n "${ROOT_DIR}/downstream_analysis/config/calibrations_template.tsv" "${DATE_DIR}/calibrations_to_finalize.tsv"

if [ ! -s "${DATE_DIR}/species_tree_with_calibrations.tre" ]; then
    echo "[ERROR] Missing finalized calibrated tree: ${DATE_DIR}/species_tree_with_calibrations.tre" >&2
    echo "[ERROR] Finalize calibrations_to_finalize.tsv and write the calibrated topology before running MCMCTree." >&2
    exit 3
fi

mkdir -p "${IQ2MC_DIR}"

set +u
conda activate "${IQTREE_ENV}"
set -u

if ! command -v iqtree3 >/dev/null 2>&1; then
    echo "[ERROR] Required command not found in IQTREE_ENV=${IQTREE_ENV}: iqtree3" >&2
    exit 2
fi

if [ ! -s "${IQ2MC_PREFIX}.mcmctree.ctl" ] || [ ! -s "${IQ2MC_PREFIX}.mcmctree.hessian" ]; then
    cd "${IQ2MC_DIR}"
    iqtree3 \
        -s "${DATE_DIR}/single_copy.4d_sites.phy" \
        -m "${IQTREE_DATING_MODEL}" \
        -te "${DATE_DIR}/species_tree_with_calibrations.tre" \
        --dating mcmctree \
        --mcmc-iter "${MCMC_BURNIN},${MCMC_SAMPFREQ},${MCMC_NSAMPLE}" \
        --mcmc-bds "${MCMC_BDS}" \
        --mcmc-clock "${MCMC_CLOCK}" \
        -T "${THREADS}" \
        --prefix "${IQ2MC_PREFIX_NAME}"
fi

set +u
conda activate "${PAML_ENV}"
set -u

if ! command -v mcmctree >/dev/null 2>&1; then
    echo "[ERROR] Required command not found in PAML_ENV=${PAML_ENV}: mcmctree" >&2
    exit 2
fi

mkdir -p "${CHAIN_DIR}"

cp "${IQ2MC_PREFIX}.mcmctree.ctl" "${CHAIN_DIR}/mcmctree.ctl"
cp "${IQ2MC_PREFIX}.mcmctree.hessian" "${CHAIN_DIR}/"
cp "${IQ2MC_PREFIX}.rooted.nwk" "${CHAIN_DIR}/"
dummy_found="no"
for dummy in "${IQ2MC_PREFIX}".dummy*; do
    if [ -s "${dummy}" ]; then
        cp "${dummy}" "${CHAIN_DIR}/"
        dummy_found="yes"
    fi
done
if [ "${dummy_found}" != "yes" ]; then
    echo "[ERROR] IQ-TREE did not generate a dummy alignment for MCMCTree under ${IQ2MC_DIR}" >&2
    exit 3
fi

cd "${CHAIN_DIR}"
mcmctree mcmctree.ctl
