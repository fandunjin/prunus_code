#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_species14_plant_outgroup_fossil55_35_r1}"
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
IQTREE_MODEL="${IQTREE_MODEL:-LG+F+G4}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -d "${TRIMMED_DIR}" ]] || die "Missing trimmed protein alignments: ${TRIMMED_DIR}"
[[ -s "${PARTITION_FILE}" ]] || die "Missing partition file: ${PARTITION_FILE}"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"
for exe in "${PYTHON_BIN}" "${IQTREE_BIN}" "${MCMCTREE_BIN}"; do
  command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done
mkdir -p "${PREP_DIR}" "${BRIDGE_DIR}" "${CHAIN_DIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_mcmctree_protein_species14.py" \
  --trimmed-dir "${TRIMMED_DIR}" --partition-file "${PARTITION_FILE}" --out-dir "${PREP_DIR}"

# Time unit is 100 Myr. These fossils constrain lineage crowns, never WGD
# or homeolog-separation ages. The 120-Ma root bound is a broad plant prior.
cat > "${OUT_ROOT}/species_tree_with_calibrations.tre" <<'EOF'
(rch,(gtr,((ppa,(psp,pzi))'>0.35<0.55',((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe))))'>0.55<1.20'));
EOF

cat > "${OUT_ROOT}/calibration_manifest.tsv" <<'EOF'
node_label	node_definition	age_ma	bound	fossil_basis	interpretation
Racemose_crown	(ppa,psp,pzi)	35	>0.35<0.55	Prunus hirsutipetala	soft minimum/upper bound for sampled racemose lineage; not a WGD age
Prunus_crown	all sampled Prunus taxa	55	>0.55<1.20	Prunus wutuensis	soft minimum/upper bound for genus-level crown proxy; not a WGD age
Root_safeguard	Rosa-Gillenia-Prunus root	120	<1.20	broad Rosaceae prior	deep-root safeguard, not a fossil calibration
EOF

cd "${BRIDGE_DIR}"
"${IQTREE_BIN}" -s "${PREP_DIR}/single_copy_protein_species14.phy" -m "${IQTREE_MODEL}" \
  -te "${OUT_ROOT}/species_tree_with_calibrations.tre" --dating mcmctree \
  --mcmc-iter "${MCMC_BURNIN},${MCMC_SAMPFREQ},${MCMC_NSAMPLE}" \
  --mcmc-bds 1,1,0.5 --mcmc-clock IND -T "${THREADS}" --prefix species14_fossil55_racemose35

ctl="${BRIDGE_DIR}/species14_fossil55_racemose35.mcmctree.ctl"
hessian="${BRIDGE_DIR}/species14_fossil55_racemose35.mcmctree.hessian"
rooted="${BRIDGE_DIR}/species14_fossil55_racemose35.rooted.nwk"
[[ -s "${ctl}" && -s "${hessian}" && -s "${rooted}" ]] || die "IQ-TREE did not produce MCMCtree bridge files"
{
  printf '14 1\n'
  printf "(rch,(gtr,((ppa,(psp,pzi))'>0.35<0.55',((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe))))'>0.55<1.20'));\n"
} > "${rooted}"
cp "${ctl}" "${CHAIN_DIR}/mcmctree.ctl"
cp "${hessian}" "${CHAIN_DIR}/in.BV"
cp "${rooted}" "${CHAIN_DIR}/"
dummy_found=no
for dummy in "${BRIDGE_DIR}"/species14_fossil55_racemose35.dummy*; do
  if [[ -s "${dummy}" ]]; then
    cp "${dummy}" "${CHAIN_DIR}/"
    dummy_found=yes
  fi
done
[[ "${dummy_found}" == yes ]] || die "MCMCtree dummy alignment missing"

sed -i \
  -e '/^[[:space:]]*hessianfile[[:space:]]*=/d' \
  -e '/^[[:space:]]*ckpfile[[:space:]]*=/d' \
  -e 's/^seqtype[[:space:]]*=.*/seqtype = 2/' \
  -e 's/^checkpoint[[:space:]]*=.*/checkpoint = 0/' \
  -e 's/^RootAge[[:space:]]*=.*/RootAge = <1.2/' \
  -e 's/^BDparas[[:space:]]*=.*/BDparas = 1 1 0.5 C/' \
  "${CHAIN_DIR}/mcmctree.ctl"

cat > "${CHAIN_DIR}/RUN_INFO.txt" <<EOF
14-tip species-level companion tree with exact requested Prunus sampling
tips: ppa, psp, pzi, pmal, pav, pspe, pcam, psal, par, pmu, pdu, ppe, gtr, rch
outgroups: rch outer; gtr sister to Prunus; no animal taxa
polyploid policy: A representative only for ppa/psp/pzi; B copies retained in the 17-tip supplementary tree
calibrations: P. hirsutipetala 35 Ma on racemose crown; P. wutuensis 55 Ma on sampled Prunus crown proxy
root safeguard: RootAge <1.2 (120 Ma), broad plant prior, not a fossil calibration
sequence type: amino acid; model=${IQTREE_MODEL}; existing trimmed single-copy families
burnin=${MCMC_BURNIN} sampfreq=${MCMC_SAMPFREQ} nsample=${MCMC_NSAMPLE}
EOF

cd "${CHAIN_DIR}"
"${MCMCTREE_BIN}" mcmctree.ctl
printf '[OK] 14-tip species-level plant-outgroup MCMCtree run completed: %s\n' "${OUT_ROOT}"
