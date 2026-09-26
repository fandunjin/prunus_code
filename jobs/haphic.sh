#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
conda activate haphic

haphic_dir="${HOME}/cx/gitsoftware/HapHiC"

fcs_clean_dir="${PROJECT_ROOT}/fcs/final_clean"
hic_data_dir="${PROJECT_ROOT}/HIC"

output_base="${PROJECT_ROOT}/haphic"

p_ctg_fa="ONT.asm.HD-14-2DNA.hic.p_ctg.final_clean.fasta"
hap1_fa="ONT.asm.HD-14-2DNA.hic.hap1.p_ctg.final_clean.fasta"
hap2_fa="ONT.asm.HD-14-2DNA.hic.hap2.p_ctg.final_clean.fasta"

# Hi-C R1/R2
hic_r1="${hic_data_dir}/H260100013A_1.fq.gz"
hic_r2="${hic_data_dir}/H260100013A_2.fq.gz"

nchrs="16"
nchrs2="32"


mkdir -p ${output_base}/p_ctg
mkdir -p ${output_base}/merge

echo "Starting p_ctg scaffolding (n=16)..."
cd ${output_base}/p_ctg

ln -sf ${fcs_clean_dir}/${p_ctg_fa} .

bwa index ${p_ctg_fa}
bwa mem -5SP -t 48 ${p_ctg_fa} ${hic_r1} ${hic_r2} | samblaster | samtools view - -@ 12 -S -h -b -F 3340 -o HiC.bam
"${haphic_dir}/utils/filter_bam" HiC.bam 1 --nm 3 --threads 36 | samtools view - -b -@ 12 -o HiC.filtered.bam

"${haphic_dir}/haphic" pipeline ${p_ctg_fa} HiC.filtered.bam ${nchrs} --threads 48 --processes 48


echo "Starting merge scaffolding (2n=32)..."
cd ${output_base}/merge

cat ${fcs_clean_dir}/${hap1_fa} ${fcs_clean_dir}/${hap2_fa} > merge.fasta

bwa index merge.fasta
bwa mem -5SP -t 48 merge.fasta ${hic_r1} ${hic_r2} | samblaster | samtools view - -@ 12 -S -h -b -F 3340 -o HiC.bam
"${haphic_dir}/utils/filter_bam" HiC.bam 1 --nm 3 --threads 36 | samtools view - -b -@ 12 -o HiC.filtered.bam

"${haphic_dir}/haphic" pipeline merge.fasta HiC.filtered.bam ${nchrs2} --threads 48 --processes 48

echo "HapHiC Pipeline Finished. Results are in ${output_base}"