#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_fossil55_racemose35_prunus_only}"
TRIMMED_DIR="${TRIMMED_DIR:-${ROOT_DIR}/downstream_analysis/phylogenetics/04_alignments/main/trimmed}"
PARTITION_FILE="${PARTITION_FILE:-${ROOT_DIR}/downstream_analysis/phylogenetics/04_alignments/main/qc/concat_partitions.nex}"
PREP_DIR="${OUT_ROOT}/01_prepared"
BRIDGE_DIR="${OUT_ROOT}/02_iqtree_mcmctree"
CHAIN_DIR="${OUT_ROOT}/03_mcmctree_chain1"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-iqtree3}"
MCMCTREE_BIN="${MCMCTREE_BIN:-mcmctree}"
THREADS="${THREADS:-8}"
MCMC_BURNIN="${MCMC_BURNIN:-200000}"
MCMC_SAMPFREQ="${MCMC_SAMPFREQ:-50}"
MCMC_NSAMPLE="${MCMC_NSAMPLE:-50000}"
IQTREE_MODEL="${IQTREE_MODEL:-LG}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -d "${TRIMMED_DIR}" ]] || die "Missing trimmed protein alignments: ${TRIMMED_DIR}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"
for exe in "${PYTHON_BIN}" "${IQTREE_BIN}" "${MCMCTREE_BIN}"; do
  command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done
mkdir -p "${PREP_DIR}" "${BRIDGE_DIR}" "${CHAIN_DIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_mcmctree_protein_fossil55_35.py" \
  --trimmed-dir "${TRIMMED_DIR}" --partition-file "${PARTITION_FILE}" --out-dir "${PREP_DIR}"

# MCMCtree time unit: 100 Myr = 1. The fossil ages are soft minimum bounds.
# P. wutuensis (55 Ma) constrains the Prunus crown; P. hirsutipetala (35 Ma)
# constrains the racemose crown (ppa,(psp,pzi)). These are not WGD ages.
cat > "${OUT_ROOT}/species_tree_with_calibrations.tre" <<'EOF'
((ppa,(psp,pzi))'>0.35<0.55',((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe))))'>0.55<1.50';
EOF

cat > "${OUT_ROOT}/calibration_manifest.tsv" <<'EOF'
node_label	node_definition	age_ma	bound	fossil_basis	interpretation
Prunus_crown	all Prunus tips (dating-tree root)	55	>0.55	Prunus wutuensis	minimum species-divergence age; not WGD age
Racemose_crown	ppa,(psp,pzi)	35	>0.35	Prunus hirsutipetala	minimum species-divergence age; not WGD age
EOF

cd "${BRIDGE_DIR}"
"${IQTREE_BIN}" -s "${PREP_DIR}/single_copy_protein.phy" -m "${IQTREE_MODEL}" \
  -te "${OUT_ROOT}/species_tree_with_calibrations.tre" --dating mcmctree \
  --mcmc-iter "${MCMC_BURNIN},${MCMC_SAMPFREQ},${MCMC_NSAMPLE}" \
  --mcmc-bds 1,1,0.5 --mcmc-clock IND -T "${THREADS}" --prefix fossil55_racemose35

ctl="${BRIDGE_DIR}/fossil55_racemose35.mcmctree.ctl"
hessian="${BRIDGE_DIR}/fossil55_racemose35.mcmctree.hessian"
rooted="${BRIDGE_DIR}/fossil55_racemose35.rooted.nwk"
[[ -s "${ctl}" && -s "${hessian}" ]] || die "IQ-TREE did not produce MCMCtree bridge files"
[[ -s "${rooted}" ]] || die "IQ-TREE did not produce rooted tree"
cp "${ctl}" "${CHAIN_DIR}/mcmctree.ctl"
cp "${hessian}" "${CHAIN_DIR}/"
cp "${rooted}" "${CHAIN_DIR}/"
dummy_found=no
for dummy in "${BRIDGE_DIR}"/fossil55_racemose35.dummy*; do
  if [[ -s "${dummy}" ]]; then cp "${dummy}" "${CHAIN_DIR}/"; dummy_found=yes; fi
done
[[ "${dummy_found}" == yes ]] || die "MCMCtree dummy alignment missing"

# Standard PAML 4.10.10 expects the approximate-likelihood matrix as in.BV.
# IQ-TREE's hessianfile/ckpfile options are only supported by modified PAML.
cp "${hessian}" "${CHAIN_DIR}/in.BV"
sed -i \
  -e '/^[[:space:]]*hessianfile[[:space:]]*=/d' \
  -e '/^[[:space:]]*ckpfile[[:space:]]*=/d' \
  -e 's/^seqtype[[:space:]]*=.*/seqtype = 2/' \
  -e 's/^checkpoint[[:space:]]*=.*/checkpoint = 0/' \
  -e 's/^BDparas[[:space:]]*=.*/BDparas = 1 1 0.5 C/' \
  "${CHAIN_DIR}/mcmctree.ctl"

cat > "${CHAIN_DIR}/RUN_INFO.txt" <<EOF
species-level protein MCMCtree run
calibrations: P. wutuensis 55 Ma (>0.55) on Prunus crown; P. hirsutipetala 35 Ma (>0.35) on racemose crown
sequence type: amino acid (existing trimmed single-copy protein alignments)
taxa: 12 Prunus species; A representatives only for PZI/PSP/PPA
burnin=${MCMC_BURNIN} sampfreq=${MCMC_SAMPFREQ} nsample=${MCMC_NSAMPLE}
EOF

cd "${CHAIN_DIR}"
"${MCMCTREE_BIN}" mcmctree.ctl
printf '[OK] MCMCtree fossil55/racemose35 completed: %s\n' "${OUT_ROOT}"
