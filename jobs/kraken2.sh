#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate kraken2

cd ${HOME}/cx/Tea/HD-14-2DNA

mkdir -p kraken2

cd kraken2

kraken2 --db ${HOME}/cx/pluspfp/db_pluspfp/ \
    --threads 48 \
    --report 14-2_report.txt \
    --use-names \
    --output 14-2.out \
    ${HOME}/cx/Tea/HD-14-2DNA/hifiasm_output/ONT.asm.HD-14-2DNA.hic.hap1.p_ctg.fasta
