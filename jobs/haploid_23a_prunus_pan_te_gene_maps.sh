#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OG="${OG:-${ROOT_DIR}/downstream_analysis/phylogenetics/02_orthofinder/main/Results_Aug13/Orthogroups/Orthogroups.tsv}"
OUT="${OUT:-${ROOT_DIR}/analysis/prunus_pan_te/followup/geneanchor_maps}"
IDROOT="${IDROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/id_maps}"
GFFROOT="${GFFROOT:-${ROOT_DIR}/haploid_rename/gff}"
PSPGFF="${PSPGFF:-${ROOT_DIR}/data/psp/haploid/annotation/genes}"
ARGS=(--orthogroups "${OG}" --output-dir "${OUT}")
ARGS+=(--sample "pziA1,${GFFROOT}/pzi_hapA1.gff3,${IDROOT}/pzi_A.id_map.tsv")
ARGS+=(--sample "pziB1,${GFFROOT}/pzi_hapB1.gff3,${IDROOT}/pzi_B.id_map.tsv")
ARGS+=(--sample "pspA1,${PSPGFF}/hapA1_helixer_edta_masked/postprocess/psp_hapA1.edta_masked.helixer.genes.gff3,${IDROOT}/psp_A.id_map.tsv")
ARGS+=(--sample "pspB1,${PSPGFF}/hapB1_helixer_edta_masked/postprocess/psp_hapB1.edta_masked.helixer.genes.gff3,${IDROOT}/psp_B.id_map.tsv")
ARGS+=(--sample "ppaA,${ROOT_DIR}/data/ppa/haploid_rename/gff/ppa_hapA.gff3,${IDROOT}/ppa_A.id_map.tsv")
ARGS+=(--sample "ppaB,${ROOT_DIR}/data/ppa/haploid_rename/gff/ppa_hapB.gff3,${IDROOT}/ppa_B.id_map.tsv")
ARGS+=(--sample "pcam,${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.gff.gz,${IDROOT}/pcam.id_map.tsv")
"${PYTHON}" "${ROOT_DIR}/scripts/prunus_pan_te_build_geneanchor_maps.py" "${ARGS[@]}"
echo "=== gene-anchor mapping QC completed ==="
