#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate quast

cd ${HOME}/cx/Tea/HD-14-2DNA

mkdir -p quast_output

quast -t 12 -o quast_output hifiasm_output/*.fasta

