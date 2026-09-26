#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/transcriptome}"
ANNOTATION_DIR="${ANNOTATION_DIR:-${ROOT_DIR}/annotation}"
RAW_TRANS_DIR="${RAW_TRANS_DIR:-${ROOT_DIR}/rawdata/trans}"
GENOME_FA="${GENOME_FA:-${ROOT_DIR}/final_assembly/pzi_2n32.chrom/pzi_2n32.chromosomes.fa}"
if [[ ! -s "${GENOME_FA}" && -s "${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa" ]]; then GENOME_FA="${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa"; fi
GENOME_GFF="${GENOME_GFF:-}"
PACBIO_INPUT_DIR="${PACBIO_INPUT_DIR:-${ROOT_DIR}/rawdata/HD-14-2PB}"
[[ -d "${PACBIO_INPUT_DIR}" ]] || PACBIO_INPUT_DIR="${ROOT_DIR}/../rawdata/prunus/HD-14-2PB"
[[ -d "${RAW_TRANS_DIR}" ]] || RAW_TRANS_DIR="${ROOT_DIR}/../rawdata/prunus/trans"
PACBIO_FASTQ="${PACBIO_FASTQ:-${PACBIO_INPUT_DIR}/m84072_260707_083113_s2.fastq.gz}"
PACBIO_BAM="${PACBIO_BAM:-${PACBIO_INPUT_DIR}/m84072_260707_083113_s2.bam}"
THREADS="${THREADS:-24}"; ALIGN_THREADS="${ALIGN_THREADS:-15}"; SORT_THREADS="${SORT_THREADS:-8}"; PACBIO_ALIGN_THREADS="${PACBIO_ALIGN_THREADS:-15}"; SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-23}"
MEMORY="${MEMORY:-100G}"; RNA_ENV="${RNA_ENV:-trinity}"; CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"; RUN_TRINITY="${RUN_TRINITY:-1}"
CLEAN_DIR="${TRANS_DIR}/01_clean_data"; ALIGN_DIR="${TRANS_DIR}/02_hisat2_align"; ALIGN_STATS_DIR="${ALIGN_DIR}/stats"; STRINGTIE_DIR="${TRANS_DIR}/03_stringtie"; TRINITY_DIR="${TRANS_DIR}/04_trinity"; PACBIO_WORK_DIR="${TRANS_DIR}/05_pacbio"; INTEGRATED_DIR="${TRANS_DIR}/06_integrated"; EVIANN_EVIDENCE_DIR="${ANNOTATION_DIR}/genes/pzi_eviann/evidence"
IDX_PREFIX="${ALIGN_DIR}/pzi_genome_idx"; LIB_IDS=(14-2_LIBY260601004066 14-2_LIBY260601004067 14-2_LIBY260601004068)
die(){ echo "ERROR: $*" >&2; exit 1; }; require_file(){ [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }; require_cmd(){ command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }; activate_environment(){ source "$CONDA_SH"; conda activate "$RNA_ENV"; }; index_is_complete(){ local p="$1" s; for s in 1 2 3 4 5 6 7 8; do [[ -s "$p.$s.ht2" || -s "$p.$s.ht2l" ]] || return 1; done; }
TISSUE_LABELS="${TISSUE_LABELS:-unknown,unknown,unknown}"; IFS=',' read -r -a TISSUES <<< "$TISSUE_LABELS"; [[ "${#TISSUES[@]}" -eq 3 ]] || die "TISSUE_LABELS must have three labels"
LIBRARY_METADATA="${LIBRARY_METADATA:-${TRANS_DIR}/library_metadata.tsv}"
mkdir -p "$TRANS_DIR" "$CLEAN_DIR" "$ALIGN_DIR" "$ALIGN_STATS_DIR" "$STRINGTIE_DIR" "$TRINITY_DIR" "$PACBIO_WORK_DIR" "$INTEGRATED_DIR" "$EVIANN_EVIDENCE_DIR"
