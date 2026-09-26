#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate hifiasm

work_dir="${HOME}/cx/Tea/HD-14-2DNA"

name="HD-14-2DNA"

ont="HD-14-2DNA.merge.pass.fastq.gz"

hic1="${HOME}/cx/Tea/HD-14-2DNA/HIC/H260100013A_1.fq.gz"
hic2="${HOME}/cx/Tea/HD-14-2DNA/HIC/H260100013A_2.fq.gz"

cd ${work_dir}
mkdir -p hifiasm_output

hifiasm -t 64 --ont \
 --h1 "${hic1}" --h2 "${hic2}" \
 -o "hifiasm_output/ONT.asm.${name}" \
 ${ont}
