#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
OG="${OG:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13/Orthogroups/Orthogroups.tsv}"
FAM="${FAM:-${PAN_ROOT}/followup/curated_summary/curated_te_driver_families.tsv}"
OUT="${OUT:-${PAN_ROOT}/followup/geneanchor}"
TARGETS="${TARGETS:-${ROOT_DIR}/config/prunus_pan_te_geneanchor_targets_representatives.tsv}"
GENE_MAP_DIR="${GENE_MAP_DIR:-${PAN_ROOT}/followup/geneanchor_maps}"
[ -s "${OG}" ] && [ -s "${FAM}" ] && [ -s "${TARGETS}" ] || { echo "ERROR: missing OG/FAM/TARGETS" >&2; exit 1; }
[ -s "${GENE_MAP_DIR}/geneanchor_mapping_qc.tsv" ] || { echo "ERROR: run haploid_23a gene-map QC first: ${GENE_MAP_DIR}/geneanchor_mapping_qc.tsv" >&2; exit 1; }
awk -F'\t' 'NR>1 && $6+0 < 0.80 {print "ERROR: low gene mapping", $1, $6 > "/dev/stderr"; bad=1} END {exit bad}' "${GENE_MAP_DIR}/geneanchor_mapping_qc.tsv"
ARGS=(--orthogroups "${OG}" --families "${FAM}" --targets "${TARGETS}" --gene-map-dir "${GENE_MAP_DIR}" --output-dir "${OUT}")
GFF_ROOT="${GFF_ROOT:-${ROOT_DIR}/haploid_rename/gff}"
ARGS+=(--query "pziA1,${GFF_ROOT}/pzi_hapA1.gff3,${PAN_ROOT}/mask/pziA1/pziA1.fa.out,${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/id_maps/pzi_A.id_map.tsv,pzi_A")
ARGS+=(--query "pziB1,${GFF_ROOT}/pzi_hapB1.gff3,${PAN_ROOT}/mask/pziB1/pziB1.fa.out,${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/id_maps/pzi_B.id_map.tsv,pzi_B")
"${PYTHON}" "${ROOT_DIR}/scripts/prunus_pan_te_geneanchor.py" "${ARGS[@]}"
echo "=== gene-anchor TE analysis completed ==="
