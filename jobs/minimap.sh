#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate minimap

mkdir ../minimap/
minimap2 -x asm5 -t 36 ${PROJECT_ROOT}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna ${PROJECT_ROOT}/hifiasm_output/ONT.asm.HD-14-2DNA.hic.hap1.p_ctg.fasta > ${PROJECT_ROOT}/minimap/pzi_vs_peach.paf