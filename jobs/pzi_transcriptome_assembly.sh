#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# After this job completes: jsub < jobs/pzi_32_helixer.sh
# Then: jsub < jobs/pzi_annotation_eviann.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/transcriptome}"
ANNOTATION_DIR="${ANNOTATION_DIR:-${ROOT_DIR}/annotation}"
RAW_TRANS_DIR="${RAW_TRANS_DIR:-${ROOT_DIR}/rawdata/trans}"
GENOME_FA="${GENOME_FA:-${ROOT_DIR}/final_assembly/pzi_2n32.chrom/pzi_2n32.chromosomes.fa}"
GENOME_GFF="${GENOME_GFF:-}"

# The project has appeared with both layouts below; resolve the PacBio input
# directory automatically while allowing an explicit PACBIO_INPUT_DIR override.
PACBIO_INPUT_DIR="${PACBIO_INPUT_DIR:-}"
if [[ -z "${PACBIO_INPUT_DIR}" ]]; then
    if [[ -d "${ROOT_DIR}/rawdata/HD-14-2PB" ]]; then
        PACBIO_INPUT_DIR="${ROOT_DIR}/rawdata/HD-14-2PB"
    elif [[ -d "${ROOT_DIR}/../rawdata/prunus/HD-14-2PB" ]]; then
        PACBIO_INPUT_DIR="${ROOT_DIR}/../rawdata/prunus/HD-14-2PB"
    else
        PACBIO_INPUT_DIR="${ROOT_DIR}/rawdata/HD-14-2PB"
    fi
fi
if [[ ! -d "${RAW_TRANS_DIR}" && -d "${ROOT_DIR}/../rawdata/prunus/trans" ]]; then
    RAW_TRANS_DIR="${ROOT_DIR}/../rawdata/prunus/trans"
fi
PACBIO_BAM="${PACBIO_BAM:-${PACBIO_INPUT_DIR}/m84072_260707_083113_s2.bam}"
PACBIO_FASTQ="${PACBIO_FASTQ:-${PACBIO_INPUT_DIR}/m84072_260707_083113_s2.fastq.gz}"

# The supplied PacBio reads are confirmed full-length Iso-Seq reads. Set ignore
# only for a transcriptome-only rerun. External collapsed isoforms may be supplied
# as both FASTA and GTF; otherwise the FASTQ is mapped and modeled below.
PACBIO_MODE="${PACBIO_MODE:-isoseq}"
PACBIO_ISOFORM_FA="${PACBIO_ISOFORM_FA:-}"
PACBIO_ISOFORM_GTF="${PACBIO_ISOFORM_GTF:-}"
PACBIO_ALIGNER="${PACBIO_ALIGNER:-minimap2}"
RUN_TRINITY="${RUN_TRINITY:-1}"
TISSUE_LABELS="${TISSUE_LABELS:-unknown,unknown,unknown}"
UPDATE_METADATA="${UPDATE_METADATA:-0}"

THREADS="${THREADS:-24}"
SORT_THREADS="${SORT_THREADS:-8}"
ALIGN_THREADS="${ALIGN_THREADS:-15}"
PACBIO_ALIGN_THREADS="${PACBIO_ALIGN_THREADS:-${ALIGN_THREADS}}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-23}"
MEMORY="${MEMORY:-100G}"
RNA_ENV="${RNA_ENV:-trinity}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

CLEAN_DIR="${TRANS_DIR}/01_clean_data"
ALIGN_DIR="${TRANS_DIR}/02_hisat2_align"
ALIGN_STATS_DIR="${ALIGN_DIR}/stats"
STRINGTIE_DIR="${TRANS_DIR}/03_stringtie"
TRINITY_DIR="${TRANS_DIR}/04_trinity"
PACBIO_WORK_DIR="${TRANS_DIR}/05_pacbio"
INTEGRATED_DIR="${TRANS_DIR}/06_integrated"
EVIANN_EVIDENCE_DIR="${ANNOTATION_DIR}/genes/pzi_eviann/evidence"
LOG_DIR="${TRANS_DIR}/logs"
IDX_PREFIX="${ALIGN_DIR}/pzi_genome_idx"

LIB_IDS=(
    "14-2_LIBY260601004066"
    "14-2_LIBY260601004067"
    "14-2_LIBY260601004068"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

IFS=',' read -r -a TISSUES <<< "${TISSUE_LABELS}"
[[ "${#TISSUES[@]}" -eq "${#LIB_IDS[@]}" ]] || \
    die "TISSUE_LABELS must contain exactly ${#LIB_IDS[@]} comma-separated labels"
LIBRARY_METADATA="${LIBRARY_METADATA:-${TRANS_DIR}/library_metadata.tsv}"

index_is_complete() {
    local prefix="$1"
    local suffix
    for suffix in 1 2 3 4 5 6 7 8; do
        [[ -s "${prefix}.${suffix}.ht2" || -s "${prefix}.${suffix}.ht2l" ]] || return 1
    done
}

activate_environment() {
    if [[ -n "${RNA_ENV}" ]]; then
        [[ -r "${CONDA_SH}" ]] || die "Conda initialization script not readable: ${CONDA_SH}"
        # shellcheck disable=SC1090
        source "${CONDA_SH}"
        conda activate "${RNA_ENV}"
    fi
}

mkdir -p "${TRANS_DIR}" "${CLEAN_DIR}" "${ALIGN_DIR}" "${ALIGN_STATS_DIR}" "${STRINGTIE_DIR}" \
         "${TRINITY_DIR}" "${PACBIO_WORK_DIR}" "${INTEGRATED_DIR}" \
         "${EVIANN_EVIDENCE_DIR}" "${LOG_DIR}"

require_file "${GENOME_FA}" "PZI reference FASTA"
if [[ -n "${GENOME_GFF}" ]]; then
    require_file "${GENOME_GFF}" "optional genome annotation GFF"
fi
case "${PACBIO_MODE}" in
    isoseq|ignore) ;;
    *) die "PACBIO_MODE must be isoseq or ignore" ;;
esac
case "${RUN_TRINITY}" in
    0|1) ;;
    *) die "RUN_TRINITY must be 0 or 1" ;;
esac
case "${UPDATE_METADATA}" in
    0|1) ;;
    *) die "UPDATE_METADATA must be 0 or 1" ;;
esac
if [[ -n "${PACBIO_ISOFORM_FA}" || -n "${PACBIO_ISOFORM_GTF}" ]]; then
    [[ -n "${PACBIO_ISOFORM_FA}" && -n "${PACBIO_ISOFORM_GTF}" ]] \
        || die "PACBIO_ISOFORM_FA and PACBIO_ISOFORM_GTF must be supplied together"
    require_file "${PACBIO_ISOFORM_FA}" "processed PacBio Iso-Seq isoform FASTA"
    require_file "${PACBIO_ISOFORM_GTF}" "processed PacBio Iso-Seq isoform GTF"
fi
for lib in "${LIB_IDS[@]}"; do
    require_file "${RAW_TRANS_DIR}/${lib}/${lib}_1.fq.gz" "${lib} raw R1"
    require_file "${RAW_TRANS_DIR}/${lib}/${lib}_2.fq.gz" "${lib} raw R2"
done

if [[ ! -s "${LIBRARY_METADATA}" || "${UPDATE_METADATA}" == "1" ]]; then
    {
        printf 'library_id\ttissue\tr1\tr2\n'
        for i in "${!LIB_IDS[@]}"; do
            printf '%s\t%s\t%s\t%s\n' \
                "${LIB_IDS[$i]}" "${TISSUES[$i]}" \
                "${RAW_TRANS_DIR}/${LIB_IDS[$i]}/${LIB_IDS[$i]}_1.fq.gz" \
                "${RAW_TRANS_DIR}/${LIB_IDS[$i]}/${LIB_IDS[$i]}_2.fq.gz"
        done
    } > "${LIBRARY_METADATA}"
fi

activate_environment
for tool in fastp hisat2-build hisat2 samtools stringtie gffread; do
    require_cmd "${tool}"
done
if [[ "${RUN_TRINITY}" == "1" ]]; then
    require_cmd Trinity
fi
if [[ "${PACBIO_MODE}" == "isoseq" && -z "${PACBIO_ISOFORM_FA}" ]]; then
    require_file "${PACBIO_FASTQ}" "PacBio full-length Iso-Seq FASTQ"
    [[ "${PACBIO_ALIGNER}" == "minimap2" ]] || die "PACBIO_ALIGNER must currently be minimap2"
    require_cmd minimap2
fi
if [[ -s "${PACBIO_BAM}" ]]; then
    samtools quickcheck -v "${PACBIO_BAM}"
fi

echo "=== PZI genome-guided transcript evidence ==="
echo "Reference       : ${GENOME_FA}"
echo "Transcript out  : ${TRANS_DIR}"
echo "Annotation input: ${EVIANN_EVIDENCE_DIR}"
echo "PacBio input    : ${PACBIO_INPUT_DIR}"
echo "PacBio mode     : ${PACBIO_MODE}"
echo "Tissue labels   : ${TISSUE_LABELS}"
echo "Library sheet   : ${LIBRARY_METADATA}"

# 1. Illumina QC, preserving each library as an independent tissue/library.
for lib in "${LIB_IDS[@]}"; do
    r1="${CLEAN_DIR}/${lib}_clean_R1.fq.gz"
    r2="${CLEAN_DIR}/${lib}_clean_R2.fq.gz"
    if [[ ! -s "${r1}" || ! -s "${r2}" || ! -s "${CLEAN_DIR}/${lib}.fastp.json" ]]; then
        fastp -i "${RAW_TRANS_DIR}/${lib}/${lib}_1.fq.gz" \
              -I "${RAW_TRANS_DIR}/${lib}/${lib}_2.fq.gz" \
              -o "${r1}" -O "${r2}" \
              -h "${CLEAN_DIR}/${lib}.fastp.html" \
              -j "${CLEAN_DIR}/${lib}.fastp.json" \
              -w "${THREADS}" -q 5 -n 10
    fi
    require_file "${r1}" "${lib} cleaned R1"
    require_file "${r2}" "${lib} cleaned R2"
done

# 2. Illumina splice-aware genome mapping and per-library BAMs.
if ! index_is_complete "${IDX_PREFIX}"; then
    hisat2-build -p "${THREADS}" "${GENOME_FA}" "${IDX_PREFIX}"
fi
for lib in "${LIB_IDS[@]}"; do
    bam="${ALIGN_DIR}/${lib}.sorted.bam"
    if [[ ! -s "${bam}" ]]; then
        tmp_bam="${bam}.tmp.$$"
        splice_sites="${ALIGN_DIR}/${lib}.novel_splicesites.tsv"
        hisat2_log="${ALIGN_STATS_DIR}/${lib}.hisat2.log"
        hisat2 -p "${ALIGN_THREADS}" -x "${IDX_PREFIX}" \
               -1 "${CLEAN_DIR}/${lib}_clean_R1.fq.gz" \
               -2 "${CLEAN_DIR}/${lib}_clean_R2.fq.gz" --dta \
               --novel-splicesite-outfile "${splice_sites}" \
               2> "${hisat2_log}" \
            | samtools sort -@ "${SORT_THREADS}" -o "${tmp_bam}" -
        mv -f "${tmp_bam}" "${bam}"
    fi
    require_file "${bam}" "${lib} sorted BAM"
    [[ -s "${bam}.bai" ]] || samtools index -@ "${SORT_THREADS}" "${bam}"
    if [[ ! -s "${ALIGN_STATS_DIR}/${lib}.flagstat.txt" ]]; then
        samtools flagstat "${bam}" > "${ALIGN_STATS_DIR}/${lib}.flagstat.txt"
    fi
    if [[ ! -s "${ALIGN_STATS_DIR}/${lib}.mapped.count" ]]; then
        samtools view -c -F 4 "${bam}" > "${ALIGN_STATS_DIR}/${lib}.mapped.count"
    fi
    mapped_count="$(tr -d '[:space:]' < "${ALIGN_STATS_DIR}/${lib}.mapped.count")"
    [[ "${mapped_count}" =~ ^[0-9]+$ ]] || die "Invalid mapped-read count for ${lib}: ${mapped_count}"
    (( mapped_count > 0 )) || die "No Illumina reads mapped for ${lib}; check GENOME_FA and read layout"
done

MERGED_BAM="${ALIGN_DIR}/all_illumina.sorted.bam"
if [[ ! -s "${MERGED_BAM}" ]]; then
    tmp_merged_bam="${MERGED_BAM}.tmp.$$"
    samtools merge -f -@ "${SAMTOOLS_THREADS}" -o "${tmp_merged_bam}" \
        "${ALIGN_DIR}/${LIB_IDS[0]}.sorted.bam" \
        "${ALIGN_DIR}/${LIB_IDS[1]}.sorted.bam" \
        "${ALIGN_DIR}/${LIB_IDS[2]}.sorted.bam"
    mv -f "${tmp_merged_bam}" "${MERGED_BAM}"
fi
[[ -s "${MERGED_BAM}.bai" ]] || samtools index -@ "${SORT_THREADS}" "${MERGED_BAM}"
if [[ ! -s "${ALIGN_STATS_DIR}/all_illumina.flagstat.txt" ]]; then
    samtools flagstat "${MERGED_BAM}" > "${ALIGN_STATS_DIR}/all_illumina.flagstat.txt"
fi

# 3. Illumina genome-guided transcript models.
STRINGTIE_GTFS=()
for lib in "${LIB_IDS[@]}"; do
    gtf="${STRINGTIE_DIR}/${lib}.gtf"
    if [[ ! -s "${gtf}" ]]; then
        if [[ -n "${GENOME_GFF}" ]]; then
            stringtie "${ALIGN_DIR}/${lib}.sorted.bam" -p "${THREADS}" \
                -G "${GENOME_GFF}" -o "${gtf}" -l "${lib}"
        else
            stringtie "${ALIGN_DIR}/${lib}.sorted.bam" -p "${THREADS}" \
                -o "${gtf}" -l "${lib}"
        fi
    fi
    require_file "${gtf}" "${lib} StringTie GTF"
    STRINGTIE_GTFS+=("${gtf}")
done
ILLUMINA_GTF="${STRINGTIE_DIR}/pzi_illumina_stringtie_merged.gtf"
if [[ ! -s "${ILLUMINA_GTF}" ]]; then
    if [[ -n "${GENOME_GFF}" ]]; then
        stringtie --merge -p "${THREADS}" -G "${GENOME_GFF}" \
            -o "${ILLUMINA_GTF}" "${STRINGTIE_GTFS[@]}"
    else
        stringtie --merge -p "${THREADS}" -o "${ILLUMINA_GTF}" \
            "${STRINGTIE_GTFS[@]}"
    fi
fi
require_file "${ILLUMINA_GTF}" "merged Illumina StringTie GTF"
ILLUMINA_FA="${STRINGTIE_DIR}/pzi_illumina_stringtie_merged.fa"
if [[ ! -s "${ILLUMINA_FA}" ]]; then
    gffread "${ILLUMINA_GTF}" -g "${GENOME_FA}" -w "${ILLUMINA_FA}"
fi
require_file "${ILLUMINA_FA}" "merged Illumina transcript FASTA"

# 4. PacBio full-length reads: map with minimap2 splice mode, then reconstruct
# long-read models with StringTie -L. If an externally collapsed Iso-Seq GTF/FASTA
# is supplied, it is used directly and the raw FASTQ is not remapped.
PACBIO_HEADER="${PACBIO_WORK_DIR}/pacbio_header.sam"
PACBIO_FLAGSTAT="${PACBIO_WORK_DIR}/pacbio_flagstat.txt"
PACBIO_REPORT="${PACBIO_WORK_DIR}/pacbio_isoseq_report.txt"
PACBIO_ALIGN_BAM="${PACBIO_WORK_DIR}/pzi_isoseq.splice.sorted.bam"
PACBIO_GTF="${PACBIO_WORK_DIR}/pzi_isoseq.stringtie.gtf"
PACBIO_FA="${PACBIO_WORK_DIR}/pzi_isoseq.stringtie.fa"
PACBIO_MODEL_METHOD="excluded"
if [[ -s "${PACBIO_BAM}" ]]; then
    if [[ ! -s "${PACBIO_HEADER}" ]]; then
        samtools view -H "${PACBIO_BAM}" > "${PACBIO_HEADER}"
    fi
    if [[ ! -s "${PACBIO_FLAGSTAT}" ]]; then
        samtools flagstat "${PACBIO_BAM}" > "${PACBIO_FLAGSTAT}"
    fi
fi

if [[ "${PACBIO_MODE}" == "isoseq" ]]; then
    if [[ -n "${PACBIO_ISOFORM_FA}" ]]; then
        PACBIO_GTF="${PACBIO_ISOFORM_GTF}"
        PACBIO_FA="${PACBIO_ISOFORM_FA}"
        PACBIO_MODEL_METHOD="external_collapsed_isoforms"
    else
        if [[ ! -s "${PACBIO_ALIGN_BAM}" ]]; then
            tmp_pb_bam="${PACBIO_ALIGN_BAM}.tmp.$$"
            minimap2 -ax splice:hq -uf --secondary=no -t "${PACBIO_ALIGN_THREADS}" \
                "${GENOME_FA}" "${PACBIO_FASTQ}" \
                | samtools sort -@ "${SORT_THREADS}" -o "${tmp_pb_bam}" -
            mv -f "${tmp_pb_bam}" "${PACBIO_ALIGN_BAM}"
        fi
        require_file "${PACBIO_ALIGN_BAM}" "PacBio splice-aligned BAM"
        [[ -s "${PACBIO_ALIGN_BAM}.bai" ]] || samtools index -@ "${SORT_THREADS}" "${PACBIO_ALIGN_BAM}"
        if [[ ! -s "${PACBIO_WORK_DIR}/pzi_isoseq.flagstat.txt" ]]; then
            samtools flagstat "${PACBIO_ALIGN_BAM}" > "${PACBIO_WORK_DIR}/pzi_isoseq.flagstat.txt"
        fi
        if [[ ! -s "${PACBIO_WORK_DIR}/pzi_isoseq.mapped.count" ]]; then
            samtools view -c -F 4 "${PACBIO_ALIGN_BAM}" > "${PACBIO_WORK_DIR}/pzi_isoseq.mapped.count"
        fi
        pacbio_mapped_count="$(tr -d '[:space:]' < "${PACBIO_WORK_DIR}/pzi_isoseq.mapped.count")"
        [[ "${pacbio_mapped_count}" =~ ^[0-9]+$ ]] || die "Invalid PacBio mapped-read count: ${pacbio_mapped_count}"
        (( pacbio_mapped_count > 0 )) || die "No PacBio reads mapped; check GENOME_FA and Iso-Seq FASTQ"
        if [[ ! -s "${PACBIO_GTF}" ]]; then
            if [[ -n "${GENOME_GFF}" ]]; then
                stringtie -L "${PACBIO_ALIGN_BAM}" -p "${THREADS}" \
                    -G "${GENOME_GFF}" -o "${PACBIO_GTF}" -l pziIsoSeq
            else
                stringtie -L "${PACBIO_ALIGN_BAM}" -p "${THREADS}" \
                    -o "${PACBIO_GTF}" -l pziIsoSeq
            fi
        fi
        require_file "${PACBIO_GTF}" "PacBio StringTie long-read GTF"
        if [[ ! -s "${PACBIO_FA}" ]]; then
            gffread "${PACBIO_GTF}" -g "${GENOME_FA}" -w "${PACBIO_FA}"
        fi
        require_file "${PACBIO_FA}" "PacBio transcript FASTA"
        PACBIO_MODEL_METHOD="minimap2_splice_plus_stringtie_L"
    fi
else
    PACBIO_GTF=""
    PACBIO_FA=""
fi

{
    echo "PZI PacBio Iso-Seq processing report"
    echo "Source BAM: ${PACBIO_BAM}"
    echo "Source FASTQ: ${PACBIO_FASTQ}"
    echo "Mode: ${PACBIO_MODE}"
    echo "Model method: ${PACBIO_MODEL_METHOD}"
    echo "Splice BAM: ${PACBIO_ALIGN_BAM}"
    echo "Transcript GTF: ${PACBIO_GTF:-none}"
    echo "Transcript FASTA: ${PACBIO_FA:-none}"
    echo "Header: ${PACBIO_HEADER}"
    echo "Flagstat: ${PACBIO_FLAGSTAT}"
    echo
    echo "The reads were mapped to the assembled PZI genome with minimap2 splice mode."
    echo "StringTie -L was used for long-read transcript model reconstruction unless external collapsed Iso-Seq files were supplied."
} > "${PACBIO_REPORT}"

if [[ -n "${PACBIO_GTF}" ]]; then
    PACBIO_INPUT_GTFS=("${ILLUMINA_GTF}" "${PACBIO_GTF}")
else
    PACBIO_INPUT_GTFS=("${ILLUMINA_GTF}")
fi

# 5. Merge short- and long-read transcript structures in genome coordinates.
ALL_GTF="${INTEGRATED_DIR}/pzi_illumina_pacbio_merged.gtf"
if [[ ! -s "${ALL_GTF}" ]]; then
    if [[ -n "${GENOME_GFF}" ]]; then
        stringtie --merge -p "${THREADS}" -G "${GENOME_GFF}" \
            -o "${ALL_GTF}" "${PACBIO_INPUT_GTFS[@]}"
    else
        stringtie --merge -p "${THREADS}" -o "${ALL_GTF}" \
            "${PACBIO_INPUT_GTFS[@]}"
    fi
fi
require_file "${ALL_GTF}" "merged Illumina/PacBio transcript GTF"
ALL_FA="${INTEGRATED_DIR}/pzi_illumina_pacbio_merged.fa"
if [[ ! -s "${ALL_FA}" ]]; then
    gffread "${ALL_GTF}" -g "${GENOME_FA}" -w "${ALL_FA}"
fi
require_file "${ALL_FA}" "merged Illumina/PacBio transcript FASTA"

# 6. Trinity is a supplementary de novo transcript source, following the
# reference paper. It is appended to the genome-guided evidence only when the
# default RUN_TRINITY=1 (set RUN_TRINITY=0 to skip this resource-intensive step).
TRINITY_FASTA="${TRINITY_DIR}/trinity_out/Trinity.fasta"
if [[ "${RUN_TRINITY}" == "1" ]]; then
    LEFT_FILES="${CLEAN_DIR}/${LIB_IDS[0]}_clean_R1.fq.gz,${CLEAN_DIR}/${LIB_IDS[1]}_clean_R1.fq.gz,${CLEAN_DIR}/${LIB_IDS[2]}_clean_R1.fq.gz"
    RIGHT_FILES="${CLEAN_DIR}/${LIB_IDS[0]}_clean_R2.fq.gz,${CLEAN_DIR}/${LIB_IDS[1]}_clean_R2.fq.gz,${CLEAN_DIR}/${LIB_IDS[2]}_clean_R2.fq.gz"
    TRINITY_OUT="${TRINITY_DIR}/trinity_out"
    if [[ ! -s "${TRINITY_FASTA}" ]]; then
        existing_trinity_entry=""
        if [[ -d "${TRINITY_OUT}" ]]; then
            existing_trinity_entry="$(find "${TRINITY_OUT}" -mindepth 1 -maxdepth 1 -print -quit)"
        fi
        [[ -z "${existing_trinity_entry}" ]] \
            || die "Trinity output exists but ${TRINITY_FASTA} is missing"
        Trinity --seqType fq --left "${LEFT_FILES}" --right "${RIGHT_FILES}" \
            --CPU "${THREADS}" --max_memory "${MEMORY}" \
            --output "${TRINITY_OUT}"
    fi
    require_file "${TRINITY_FASTA}" "optional Trinity FASTA"
else
    TRINITY_FASTA="not_run"
fi

# 7. Stage transcript evidence for the EviAnn job under annotation/. Rebuild the
# staged FASTA each run so changing RUN_TRINITY cannot leave stale evidence.
EVIANN_STAGED_FA="${EVIANN_EVIDENCE_DIR}/pzi_illumina_pacbio.transcripts.fa"
EVIANN_STAGED_FA_TMP="${EVIANN_STAGED_FA}.tmp.$$"
cp -f "${ALL_FA}" "${EVIANN_STAGED_FA_TMP}"
if [[ "${RUN_TRINITY}" == "1" ]]; then
    cat "${TRINITY_FASTA}" >> "${EVIANN_STAGED_FA_TMP}"
fi
mv -f "${EVIANN_STAGED_FA_TMP}" "${EVIANN_STAGED_FA}"
cp -f "${ALL_GTF}" "${EVIANN_EVIDENCE_DIR}/pzi_illumina_pacbio.transcripts.gtf"
if [[ "${RUN_TRINITY}" == "1" ]]; then
    cp -f "${TRINITY_FASTA}" "${EVIANN_EVIDENCE_DIR}/trinity.supplementary.transcripts.fa"
fi
for lib in "${LIB_IDS[@]}"; do
    ln -sfnT "${ALIGN_DIR}/${lib}.sorted.bam" "${EVIANN_EVIDENCE_DIR}/${lib}.sorted.bam"
    [[ -s "${ALIGN_DIR}/${lib}.sorted.bam.bai" ]] && \
        ln -sfnT "${ALIGN_DIR}/${lib}.sorted.bam.bai" "${EVIANN_EVIDENCE_DIR}/${lib}.sorted.bam.bai"
done
ln -sfnT "${MERGED_BAM}" "${EVIANN_EVIDENCE_DIR}/all_illumina.sorted.bam"
[[ -s "${MERGED_BAM}.bai" ]] && \
    ln -sfnT "${MERGED_BAM}.bai" "${EVIANN_EVIDENCE_DIR}/all_illumina.sorted.bam.bai"
if [[ -s "${PACBIO_ALIGN_BAM}" ]]; then
    ln -sfnT "${PACBIO_ALIGN_BAM}" "${EVIANN_EVIDENCE_DIR}/pzi_isoseq.splice.sorted.bam"
    [[ -s "${PACBIO_ALIGN_BAM}.bai" ]] && ln -sfnT "${PACBIO_ALIGN_BAM}.bai" "${EVIANN_EVIDENCE_DIR}/pzi_isoseq.splice.sorted.bam.bai"
fi

{
    echo "PZI genome-guided transcript evidence manifest"
    echo "Generated: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Reference FASTA: ${GENOME_FA}"
    echo "Raw Illumina directory: ${RAW_TRANS_DIR}"
    printf 'Illumina libraries: %s\n' "${LIB_IDS[*]}"
    echo "Library metadata: ${LIBRARY_METADATA}"
    echo "Tissue labels: ${TISSUE_LABELS}"
    echo "Merged Illumina BAM: ${MERGED_BAM}"
    echo "Illumina StringTie GTF: ${ILLUMINA_GTF}"
    echo "PacBio source BAM: ${PACBIO_BAM}"
    echo "PacBio source FASTQ: ${PACBIO_FASTQ}"
    echo "PacBio model method: ${PACBIO_MODEL_METHOD}"
    echo "PacBio alignment threads: ${PACBIO_ALIGN_THREADS}"
    echo "PacBio processing report: ${PACBIO_REPORT}"
    echo "PacBio splice BAM: ${PACBIO_ALIGN_BAM}"
    echo "PacBio transcript GTF: ${PACBIO_GTF:-none}"
    echo "Merged transcript GTF: ${ALL_GTF}"
    echo "Merged transcript FASTA: ${ALL_FA}"
    echo "Trinity FASTA: ${TRINITY_FASTA}"
    echo "Trinity included in EviAnn transcript evidence: ${RUN_TRINITY}"
    echo "EviAnn evidence directory: ${EVIANN_EVIDENCE_DIR}"
    echo
    echo "The merged GTF/FASTA are genome-guided transcript evidence, not final gene annotation."
    echo "Run Helixer first and then jobs/pzi_annotation_eviann.sh for the annotation GFF3/CDS/protein."
} > "${TRANS_DIR}/transcriptome_manifest.txt"

echo "=== Transcript evidence completed ==="
echo "Merged GTF  : ${ALL_GTF}"
echo "Merged FASTA: ${ALL_FA}"
echo "EviAnn input: ${EVIANN_EVIDENCE_DIR}/pzi_illumina_pacbio.transcripts.fa"
