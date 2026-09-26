#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
TE_GFF3="${TE_GFF3:-${ROOT_DIR}/annotation/chr32_edta/pzi_2n32.chromosomes.fa.mod.EDTA.TEanno.gff3}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/annotation/genes/pzi_eviann/postprocess/structure_qc}"
BAM_SPECS=(
  "illumina_066=${ROOT_DIR}/transcriptome/02_hisat2_align/14-2_LIBY260601004066.sorted.bam"
  "illumina_067=${ROOT_DIR}/transcriptome/02_hisat2_align/14-2_LIBY260601004067.sorted.bam"
  "illumina_068=${ROOT_DIR}/transcriptome/02_hisat2_align/14-2_LIBY260601004068.sorted.bam"
  "pacbio_isoseq=${ROOT_DIR}/transcriptome/05_pacbio/pzi_isoseq.splice.sorted.bam"
)
source "${ROOT_DIR}/jobs/annotation_structure_qc_common.sh"
run_structure_qc
