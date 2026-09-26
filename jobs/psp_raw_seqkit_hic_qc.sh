#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT="${ROOT}/psp_raw_qc"
HIFI_FASTQ="${ROOT}/CRR2018654.fastq.gz"
HIFI_BAM="${ROOT}/CRR2018654.bam"
HIC_R1="${ROOT}/CRR2018655_r1.fq.gz"
HIC_R2="${ROOT}/CRR2018655_r2.fq.gz"
ILL_R1="${ROOT}/CRR2018656_r1.fq.gz"
ILL_R2="${ROOT}/CRR2018656_r2.fq.gz"
MERGE_FILTERED_BAM="${ROOT}/haphic/merge/HiC.filtered.bam"
FINAL_FA="${ROOT}/haphic/merge/04.build/final_assembly_qc_outputs/psp_2n32.chromosomes.fa"
HAPHIC_DIR="${HOME}/cx/gitsoftware/HapHiC"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"

# Set REMAP_FINAL=1 to reproduce the PZI final-assembly Hi-C filtering definition.
REMAP_FINAL="${REMAP_FINAL:-0}"

source "${CONDA_SH}"
conda activate seqkit
mkdir -p "${OUT}"

for f in "${HIC_R1}" "${HIC_R2}" "${ILL_R1}" "${ILL_R2}"; do
    [[ -s "${f}" ]] || { echo "Missing FASTQ: ${f}" >&2; exit 1; }
done

if [[ -s "${HIFI_FASTQ}" ]]; then
    HIFI_INPUT="${HIFI_FASTQ}"
elif [[ -s "${HIFI_BAM}" ]]; then
    echo "CRR2018654.fastq.gz is absent; converting the public BAM to a temporary FASTQ stream."
    samtools fastq -@ 16 "${HIFI_BAM}" | gzip > "${OUT}/CRR2018654.from_bam.fastq.gz"
    HIFI_INPUT="${OUT}/CRR2018654.from_bam.fastq.gz"
else
    echo "Missing PacBio HiFi BAM/FASTQ" >&2
    exit 1
fi

seqkit stats -a -T -j 8 \
    "${HIFI_INPUT}" "${HIC_R1}" "${HIC_R2}" "${ILL_R1}" "${ILL_R2}" \
    | tee "${OUT}/psp_raw_seqkit.tsv"

awk -F '\t' '
    NR == 1 {
        print "file\tread_count\ttotal_bases\tread_N50\tQ20_percent\tQ30_percent\tQ1"
        next
    }
    {
        gsub(",", "", $4); gsub(",", "", $5)
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $4, $5, $13, $15, $16, $9
    }
' "${OUT}/psp_raw_seqkit.tsv" > "${OUT}/psp_raw_read_metrics.tsv"

echo "=== PSP raw read metrics ==="
cat "${OUT}/psp_raw_read_metrics.tsv"

if [[ -s "${MERGE_FILTERED_BAM}" ]]; then
    samtools quickcheck -v "${MERGE_FILTERED_BAM}"
    samtools flagstat -@ 16 "${MERGE_FILTERED_BAM}" > "${OUT}/merge_HiC.filtered.flagstat.txt"
    merge_pairs=$(samtools view -@ 16 -c -f 65 -F 3340 "${MERGE_FILTERED_BAM}")
    printf "merge_filtered_HiC_read_pairs\t%s\n" "${merge_pairs}" | tee "${OUT}/psp_filtered_HiC_pairs.tsv"
else
    echo "Missing existing HapHiC filtered BAM: ${MERGE_FILTERED_BAM}" >&2
fi

if [[ "${REMAP_FINAL}" == "1" ]]; then
    [[ -s "${FINAL_FA}" ]] || { echo "Missing final chromosome FASTA: ${FINAL_FA}" >&2; exit 1; }
    [[ -s "${FINAL_FA}.fai" ]] || samtools faidx "${FINAL_FA}"
    bwa index "${FINAL_FA}" 2>/dev/null || true
    raw_bam="${OUT}/psp.final32.HiC.bam"
    filtered_bam="${OUT}/psp.final32.HiC.filtered.bam"
    if [[ ! -s "${raw_bam}" ]]; then
        bwa mem -5SP -t 48 "${FINAL_FA}" "${HIC_R1}" "${HIC_R2}" \
            | samblaster \
            | samtools view -@ 16 -S -h -b -F 3340 -o "${raw_bam}"
    fi
    if [[ ! -s "${filtered_bam}" ]]; then
        "${HAPHIC_DIR}/utils/filter_bam" "${raw_bam}" 1 --nm 3 --threads 36 \
            | samtools view -@ 16 -b -o "${filtered_bam}"
    fi
    samtools flagstat -@ 16 "${filtered_bam}" > "${OUT}/psp.final32.HiC.filtered.flagstat.txt"
    chrom_pairs=$(awk 'NR == FNR { chrom[$1] = 1; next } { flag=$2; mate=$7; if ((flag & 64) && $5 >= 10 && ($3 in chrom) && (mate == "=" || mate in chrom)) n++ } END { printf "%.0f", n+0 }' \
        "${FINAL_FA}.fai" <(samtools view -@ 16 "${filtered_bam}"))
    printf "final32_filtered_HiC_read_pairs_mapq10\t%s\n" "${chrom_pairs}" | tee -a "${OUT}/psp_filtered_HiC_pairs.tsv"
fi
