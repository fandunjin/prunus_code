#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/psp}"
TRANS_DIR="${TRANS_DIR:-${DATA_DIR}/transcriptome_32}"
GENOME_FA="${GENOME_FA:-${DATA_DIR}/haphic/merge/04.build/final_assembly_qc_outputs/psp_2n32.chromosomes.fa}"
RAW_DIR="${RAW_DIR:-${DATA_DIR}}"
THREADS="${THREADS:-24}"; ALIGN_THREADS="${ALIGN_THREADS:-15}"; SORT_THREADS="${SORT_THREADS:-8}"; SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-23}"; MEMORY="${MEMORY:-100G}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"; RNA_ENV="${RNA_ENV:-trinity}"; RUN_TRINITY="${RUN_TRINITY:-1}"
CLEAN_DIR="${TRANS_DIR}/01.clean_data"; ALIGN_DIR="${TRANS_DIR}/02.hisat2_align"; TRINITY_DIR="${TRANS_DIR}/03.trinity_assembly"; STRINGTIE_DIR="${TRANS_DIR}/04.stringtie"; INTEGRATED_DIR="${TRANS_DIR}/05.integrated"
EVIANN_EVIDENCE_DIR="${DATA_DIR}/annotation/genes/psp32_eviann/evidence"; LIB_IDS=(leaf root stem); R1S=(CRR2018657_r1.fq.gz CRR2018658_r1.fq.gz CRR2018659_r1.fq.gz); R2S=(CRR2018657_r2.fq.gz CRR2018658_r2.fq.gz CRR2018659_r2.fq.gz)
IDX_PREFIX="${ALIGN_DIR}/psp_genome_idx"; die(){ echo "ERROR: $*" >&2; exit 1; }; require_file(){ [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }; require_cmd(){ command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }; activate_environment(){ source "$CONDA_SH"; conda activate "$RNA_ENV"; }; index_is_complete(){ local p="$1" s; for s in 1 2 3 4 5 6 7 8; do [[ -s "$p.$s.ht2" || -s "$p.$s.ht2l" ]] || return 1; done; }
mkdir -p "$TRANS_DIR" "$CLEAN_DIR" "$ALIGN_DIR" "$TRINITY_DIR" "$STRINGTIE_DIR" "$INTEGRATED_DIR" "$EVIANN_EVIDENCE_DIR"
