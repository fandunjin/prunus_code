#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
TE_GFF3="${TE_GFF3:-${ROOT_DIR}/data/psp/annotation/repeats/chr32_edta/psp_2n32.chromosomes.fa.mod.EDTA.TEanno.gff3}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/data/psp/annotation/genes/psp32_eviann/postprocess/structure_qc}"
BAM_SPECS=(
  "leaf=${ROOT_DIR}/data/psp/transcriptome_32/02.hisat2_align/leaf_sorted.bam"
  "root=${ROOT_DIR}/data/psp/transcriptome_32/02.hisat2_align/root_sorted.bam"
  "stem=${ROOT_DIR}/data/psp/transcriptome_32/02.hisat2_align/stem_sorted.bam"
)
source "${ROOT_DIR}/jobs/annotation_structure_qc_common.sh"
run_structure_qc
