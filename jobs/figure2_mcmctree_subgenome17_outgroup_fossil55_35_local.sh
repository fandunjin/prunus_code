#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_subgenome17_outgroup_fossil55_r1}"
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

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_mcmctree_protein_subgenome17.py" \
  --trimmed-dir "${TRIMMED_DIR}" --partition-file "${PARTITION_FILE}" --out-dir "${PREP_DIR}"

cat > "${OUT_ROOT}/species_tree_with_calibrations.tre" <<'EOF'
(rch,(gtr,((ppa_A,(psp_A,pzi_A)),((ppa_B,(psp_B,pzi_B)),((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe)))))'>0.55<1.20'));
EOF

cat > "${OUT_ROOT}/calibration_manifest.tsv" <<'EOF'
node_label	node_definition	age_ma	bound	fossil_basis	interpretation
Racemose_fossil	no monophyletic six-tip node under the ASTRAL topology	35	not applied	Prunus hirsutipetala	reserved for the species-level 14-tip tree; not forced onto homeolog tips
Prunus_crown	all sampled Prunus subgenome tips	55	>0.55<1.20	Prunus wutuensis	soft minimum/upper bound for genus-level crown proxy; not a WGD age
Root_safeguard	Rosa-Gillenia-Prunus root	120	<1.20	broad Rosaceae prior	deep-root safeguard, not a fossil calibration
EOF

cd "${BRIDGE_DIR}"
"${IQTREE_BIN}" -s "${PREP_DIR}/single_copy_protein_subgenome17.phy" -m "${IQTREE_MODEL}" \
  -te "${OUT_ROOT}/species_tree_with_calibrations.tre" --dating mcmctree \
  --mcmc-iter "${MCMC_BURNIN},${MCMC_SAMPFREQ},${MCMC_NSAMPLE}" \
  --mcmc-bds 1,1,0.5 --mcmc-clock IND -T "${THREADS}" --prefix subgenome17_fossil55_racemose35

ctl="${BRIDGE_DIR}/subgenome17_fossil55_racemose35.mcmctree.ctl"
hessian="${BRIDGE_DIR}/subgenome17_fossil55_racemose35.mcmctree.hessian"
rooted="${BRIDGE_DIR}/subgenome17_fossil55_racemose35.rooted.nwk"
[[ -s "${ctl}" && -s "${hessian}" && -s "${rooted}" ]] || die "IQ-TREE did not produce MCMCtree bridge files"
# IQ-TREE may relocate an internal calibration label when serializing a
# rooted tree.  The Hessian is independent of fossil labels, so restore the
# intended explicit node annotations before passing the tree to MCMCtree.
{
  printf '17 1\n'
  printf "(rch,(gtr,((ppa_A,(psp_A,pzi_A)),((ppa_B,(psp_B,pzi_B)),((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe)))))'>0.55<1.20'));\n"
} > "${rooted}"
cp "${ctl}" "${CHAIN_DIR}/mcmctree.ctl"
cp "${hessian}" "${CHAIN_DIR}/in.BV"
cp "${rooted}" "${CHAIN_DIR}/"
dummy_found=no
for dummy in "${BRIDGE_DIR}"/subgenome17_fossil55_racemose35.dummy*; do
  if [[ -s "${dummy}" ]]; then cp "${dummy}" "${CHAIN_DIR}/"; dummy_found=yes; fi
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
17-tip subgenome-resolved protein MCMCtree run
tips: ppa_A/B, psp_A/B, pzi_A/B, 11 diploid or plant outgroup taxa
outgroups: rch outer; gtr sister to Prunus
calibrations: P. hirsutipetala 35 Ma reserved for the 14-tip species tree; P. wutuensis 55 Ma on sampled Prunus crown proxy
root safeguard: RootAge <1.2 (120 Ma), broad plant prior, not a fossil calibration
sequence type: amino acid; model=${IQTREE_MODEL}; existing trimmed single-copy families
burnin=${MCMC_BURNIN} sampfreq=${MCMC_SAMPFREQ} nsample=${MCMC_NSAMPLE}
EOF

cd "${CHAIN_DIR}"
"${MCMCTREE_BIN}" mcmctree.ctl
printf '[OK] 17-tip subgenome/outgroup MCMCtree run completed: %s\n' "${OUT_ROOT}"
