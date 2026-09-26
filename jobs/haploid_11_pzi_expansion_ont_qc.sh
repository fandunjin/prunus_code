#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/analysis/pzi_expansion_qc}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"

# Replace with the exact raw ONT fastq used for the PZI assembly when needed.
ONT_READS="${ONT_READS:-${ROOT_DIR}/rawdata/HD-14-2DNA.pass.fastq.gz}"
# Optional directory with haplotype-binned ONT reads. When unset, the shared
# tetraploid read set is used and cross-haplotype mapping remains a caveat.
HAP_READS_DIR="${HAP_READS_DIR:-}"

# The default comparators are intended for an assembly-support screen only.
# A gzipped FASTA is accepted by minimap2. Set COMPARATORS to a space-separated
# list of label:path pairs to change the set, for example:
# PSP captures within-species haplotype diversity. PPA, PCAM and PPE retain
# representatives of the three focal Prunus inflorescence categories.
# COMPARATORS="pspA1:/.../psp_hapA1.fa pspA2:/.../psp_hapA2.fa pspB1:/.../psp_hapB1.fa pspB2:/.../psp_hapB2.fa ppa:/.../ppa.fa.gz pcam:/.../pcam.fa.gz ppe:/.../ppe.fa"
COMPARATORS="${COMPARATORS:-pspA1:${ROOT_DIR}/data/psp/haploid/genome/psp_hapA1.fa pspA2:${ROOT_DIR}/data/psp/haploid/genome/psp_hapA2.fa pspB1:${ROOT_DIR}/data/psp/haploid/genome/psp_hapB1.fa pspB2:${ROOT_DIR}/data/psp/haploid/genome/psp_hapB2.fa ppa:${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz pcam:${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz ppe:${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna}"

THREADS="${THREADS:-8}"
MAP_THREADS="${MAP_THREADS:-6}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-2}"
MIN_MAPQ="${MIN_MAPQ:-20}"
WINDOW_SIZE="${WINDOW_SIZE:-100000}"
MIN_PAF_BLOCK="${MIN_PAF_BLOCK:-5000}"
MIN_CANDIDATE_LEN="${MIN_CANDIDATE_LEN:-10000}"
MERGE_GAP="${MERGE_GAP:-2000}"
FLANK_BP="${FLANK_BP:-1000}"
MIN_SPANNING_READS="${MIN_SPANNING_READS:-3}"
MIN_DEPTH_RATIO="${MIN_DEPTH_RATIO:-0.5}"
MAX_DEPTH_RATIO="${MAX_DEPTH_RATIO:-1.5}"
FORCE_RERUN="${FORCE_RERUN:-0}"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
PLOT_ENV="${PLOT_ENV:-haphic}"
QC_SCRIPT="${QC_SCRIPT:-${ROOT_DIR}/scripts/pzi_expansion_qc.py}"

HAPS=(hapA1 hapA2 hapB1 hapB2)

# ===================== Pipeline starts here =====================

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

remove_if_forced() {
    if [ "${FORCE_RERUN}" = "1" ]; then
        rm -f "$@"
    fi
}

set_hap_reads() {
    local hap="$1"
    local candidate=""
    HAP_READS="${ONT_READS}"
    if [ -n "${HAP_READS_DIR}" ]; then
        for candidate in \
            "${HAP_READS_DIR}/pzi_${hap}.ont.fastq.gz" \
            "${HAP_READS_DIR}/pzi_${hap}.fastq.gz" \
            "${HAP_READS_DIR}/${hap}.fastq.gz"; do
            if [ -s "${candidate}" ]; then
                HAP_READS="${candidate}"
                return
            fi
        done
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${ONT_READS}" "PZI ONT reads"
require_file "${QC_SCRIPT}" "PZI expansion QC Python helper"
for hap in "${HAPS[@]}"; do
    require_file "${HAP_GENOME_DIR}/pzi_${hap}.fa" "${hap} final haplotype FASTA"
done
for label_path in ${COMPARATORS}; do
    label="${label_path%%:*}"
    fasta="${label_path#*:}"
    require_file "${fasta}" "Comparator FASTA ${label}"
done

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 was not found." >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools was not found." >&2; exit 1; }

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 8) else "ERROR: Python >= 3.8 is required.")'

echo "=== PZI genome-expansion ONT support QC ==="
echo "Haplotype FASTAs : ${HAP_GENOME_DIR}"
echo "ONT reads        : ${ONT_READS}"
echo "Comparators      : ${COMPARATORS}"
echo "Output directory : ${OUT_DIR}"
echo "Minimum MAPQ     : ${MIN_MAPQ}"
echo "Window size      : ${WINDOW_SIZE} bp"
echo

for hap in "${HAPS[@]}"; do
    conda activate "${MINIMAP_ENV}"
    genome="${HAP_GENOME_DIR}/pzi_${hap}.fa"
    hap_dir="${OUT_DIR}/${hap}"
    alignment_dir="${hap_dir}/comparative_alignment"
    depth_dir="${hap_dir}/ont_depth"
    support_dir="${hap_dir}/candidate_support"
    self_dir="${hap_dir}/self_redundancy"
    mkdir -p "${alignment_dir}" "${depth_dir}" "${support_dir}" "${self_dir}"

    if [ ! -s "${genome}.fai" ]; then
        samtools faidx "${genome}"
    fi
    fai="${genome}.fai"
    set_hap_reads "${hap}"
    reads="${HAP_READS}"
    require_file "${reads}" "ONT reads for ${hap}"

    echo "=== ${hap} ==="
    echo "Genome : ${genome}"
    echo "Reads  : ${reads}"
    if [ "${reads}" = "${ONT_READS}" ]; then
        echo "WARNING: shared tetraploid ONT reads are used; primary mapping can cross-map among homologs."
    fi

    paf_args=()
    for label_path in ${COMPARATORS}; do
        label="${label_path%%:*}"
        comparator="${label_path#*:}"
        paf="${alignment_dir}/${label}_vs_pzi_${hap}.asm20.paf"
        remove_if_forced "${paf}"
        if [ ! -s "${paf}" ]; then
            echo ">>> Whole-genome alignment: ${label} query -> PZI ${hap} target"
            minimap2 -x asm20 --secondary=no -t "${MAP_THREADS}" "${genome}" "${comparator}" > "${paf}"
        else
            echo ">>> Existing comparator PAF found: ${paf}"
        fi
        paf_args+=(--paf "${label}:${paf}")
    done

    candidates_bed="${alignment_dir}/pzi_${hap}.non_aligned_candidates.bed"
    coverage_tsv="${alignment_dir}/pzi_${hap}.comparator_alignment_coverage.tsv"
    candidate_summary="${alignment_dir}/pzi_${hap}.non_aligned_candidate_summary.tsv"
    remove_if_forced "${candidates_bed}" "${coverage_tsv}" "${candidate_summary}"
    if [ ! -s "${candidates_bed}" ]; then
        echo ">>> Calling non-aligned PZI candidate intervals"
        "${PYTHON_BIN}" "${QC_SCRIPT}" candidates \
            --fai "${fai}" \
            "${paf_args[@]}" \
            --candidates-bed "${candidates_bed}" \
            --coverage-tsv "${coverage_tsv}" \
            --candidate-summary-tsv "${candidate_summary}" \
            --min-paf-block "${MIN_PAF_BLOCK}" \
            --min-mapq "${MIN_MAPQ}" \
            --merge-gap "${MERGE_GAP}" \
            --min-candidate-len "${MIN_CANDIDATE_LEN}"
    fi

    bam="${depth_dir}/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.sorted.bam"
    bam_tmp="${bam}.tmp.bam"
    remove_if_forced "${bam}" "${bam}.bai" "${bam_tmp}"
    if [ -s "${bam}" ] && ! samtools quickcheck "${bam}"; then
        echo "WARNING: incomplete or corrupt ONT BAM will be rebuilt: ${bam}" >&2
        rm -f "${bam}" "${bam}.bai"
    fi
    if [ ! -s "${bam}" ]; then
        echo ">>> ONT read mapping: ${hap}"
        rm -f "${bam_tmp}"
        minimap2 -ax map-ont -t "${MAP_THREADS}" "${genome}" "${reads}" \
            | samtools view -@ "${SAMTOOLS_THREADS}" -b -F 2308 -q "${MIN_MAPQ}" \
            | samtools sort -@ "${SAMTOOLS_THREADS}" -o "${bam_tmp}"
        samtools quickcheck "${bam_tmp}"
        mv "${bam_tmp}" "${bam}"
    else
        echo ">>> Existing ONT BAM found: ${bam}"
    fi
    if [ ! -s "${bam}.bai" ]; then
        samtools index -@ "${SAMTOOLS_THREADS}" "${bam}"
    fi
    samtools flagstat -@ "${SAMTOOLS_THREADS}" "${bam}" > "${depth_dir}/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.flagstat.txt"

    windows_tsv="${depth_dir}/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.${WINDOW_SIZE}bp_windows.tsv"
    depth_summary="${depth_dir}/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.depth_by_chromosome.tsv"
    depth_done="${depth_dir}/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.${WINDOW_SIZE}bp_depth.complete"
    windows_tmp="${windows_tsv}.tmp"
    summary_tmp="${depth_summary}.tmp"
    remove_if_forced "${windows_tsv}" "${depth_summary}" "${depth_done}" "${windows_tmp}" "${summary_tmp}"
    if [ ! -s "${windows_tsv}" ] || [ ! -s "${depth_summary}" ] || [ ! -s "${depth_done}" ]; then
        echo ">>> Calculating ${WINDOW_SIZE}-bp ONT-depth windows: ${hap}"
        # The BAM was already restricted to primary MAPQ >= MIN_MAPQ records
        # by samtools view.  Use -Q (alignment MAPQ), not -q (base quality),
        # and avoid the unsupported samtools-depth -F option found in older
        # cluster builds.  Temporary outputs prevent a failed pipe from being
        # mistaken for a completed depth calculation on restart.
        rm -f "${windows_tmp}" "${summary_tmp}" "${depth_done}"
        samtools depth -@ "${SAMTOOLS_THREADS}" -aa -Q "${MIN_MAPQ}" "${bam}" \
            | "${PYTHON_BIN}" "${QC_SCRIPT}" depth \
                --fai "${fai}" \
                --depth-input - \
                --windows-tsv "${windows_tmp}" \
                --summary-tsv "${summary_tmp}" \
                --window-size "${WINDOW_SIZE}"
        require_file "${windows_tmp}" "ONT depth windows for ${hap}"
        require_file "${summary_tmp}" "ONT chromosome depth summary for ${hap}"
        mv "${windows_tmp}" "${windows_tsv}"
        mv "${summary_tmp}" "${depth_summary}"
        printf "completed\n" > "${depth_done}"
    fi

    support_tsv="${support_dir}/pzi_${hap}.candidate_ont_support.tsv"
    passing_tsv="${support_dir}/pzi_${hap}.candidate_ont_supported.tsv"
    support_summary="${support_dir}/pzi_${hap}.candidate_ont_support_summary.tsv"
    remove_if_forced "${support_tsv}" "${passing_tsv}" "${support_summary}"
    if [ ! -s "${support_tsv}" ]; then
        echo ">>> Counting ONT reads crossing both candidate boundaries: ${hap}"
        "${PYTHON_BIN}" "${QC_SCRIPT}" support \
            --bam "${bam}" \
            --candidates-bed "${candidates_bed}" \
            --depth-windows "${windows_tsv}" \
            --output-tsv "${support_tsv}" \
            --passing-tsv "${passing_tsv}" \
            --summary-tsv "${support_summary}" \
            --flank-bp "${FLANK_BP}" \
            --min-spanning-reads "${MIN_SPANNING_READS}" \
            --min-mapq "${MIN_MAPQ}" \
            --min-depth-ratio "${MIN_DEPTH_RATIO}" \
            --max-depth-ratio "${MAX_DEPTH_RATIO}"
    fi

    figure_prefix="${hap_dir}/figures/pzi_${hap}"
    mkdir -p "${hap_dir}/figures"
    if [ "${FORCE_RERUN}" = "1" ] || [ ! -s "${figure_prefix}.ont_depth_candidates.pdf" ]; then
        echo ">>> Drawing ONT-depth, comparative-coverage and candidate-support QC plots: ${hap}"
        conda activate "${PLOT_ENV}"
        if command -v python3 >/dev/null 2>&1; then
            PLOT_PYTHON_BIN="python3"
        else
            PLOT_PYTHON_BIN="python"
        fi
        "${PLOT_PYTHON_BIN}" -c 'import matplotlib; print("matplotlib:", matplotlib.__version__)'
        "${PLOT_PYTHON_BIN}" "${QC_SCRIPT}" visualize \
            --depth-windows "${windows_tsv}" \
            --candidates-bed "${candidates_bed}" \
            --coverage-tsv "${coverage_tsv}" \
            --support-tsv "${support_tsv}" \
            --out-prefix "${figure_prefix}" \
            --title "P. zippeliana ${hap}" \
            --min-depth-ratio "${MIN_DEPTH_RATIO}" \
            --max-depth-ratio "${MAX_DEPTH_RATIO}"
    fi

    self_paf="${self_dir}/pzi_${hap}.self.asm5.paf"
    self_summary="${self_dir}/pzi_${hap}.cross_chromosome_alignment_summary.tsv"
    remove_if_forced "${self_paf}" "${self_summary}"
    if [ ! -s "${self_paf}" ]; then
        echo ">>> Self-alignment to screen cross-chromosome redundancy: ${hap}"
        minimap2 -x asm5 -D -t "${MAP_THREADS}" "${genome}" "${genome}" > "${self_paf}"
    fi
    if [ ! -s "${self_summary}" ]; then
        awk -v min_block="${MIN_PAF_BLOCK}" -v min_mapq="${MIN_MAPQ}" 'BEGIN { OFS="\t" }
            $1 != $6 && $11 >= min_block && $12 >= min_mapq {
                key=$1 SUBSEP $6; aligned[key]+=$11; matching[key]+=$10; records[key]++
            }
            END {
                print "query_seq","target_seq","aligned_bp","matching_bp","alignment_records"
                for (key in aligned) {
                    split(key,a,SUBSEP); print a[1],a[2],aligned[key],matching[key],records[key]
                }
            }' "${self_paf}" | sort -k1,1 -k3,3nr > "${self_summary}"
    fi

    {
        printf "PZI expansion assembly-support review: %s\n" "${hap}"
        printf "================================================\n\n"
        printf "Interpretation rules:\n"
        printf -- "- Non-aligned candidates are the complement of filtered whole-genome PAF blocks from the selected comparators. They are not automatically PZI-specific insertions.\n"
        printf -- "- A candidate passes ONT support only when the same ONT reads cross both boundaries with %s-bp flanks; at least %s primary MAPQ >= %s reads are required and its window-depth ratio is %.2f-%.2f.\n" "${FLANK_BP}" "${MIN_SPANNING_READS}" "${MIN_MAPQ}" "${MIN_DEPTH_RATIO}" "${MAX_DEPTH_RATIO}"
        printf -- "- When shared tetraploid ONT reads are used, inspect representative loci with haplotype-binned reads or unique flanks before making a definitive claim.\n"
        printf -- "- Review cross-chromosome self-alignments for abnormal near-full-length duplicated blocks; do not interpret expected between-haplotype homology as redundancy in this within-haplotype screen.\n\n"
        printf "Core output files:\n"
        printf -- "- %s\n" "${candidate_summary}"
        printf -- "- %s\n" "${depth_summary}"
        printf -- "- %s\n" "${support_summary}"
        printf -- "- %s\n" "${passing_tsv}"
        printf -- "- %s\n" "${self_summary}"
        printf -- "- %s.ont_depth_candidates.pdf\n" "${figure_prefix}"
        printf -- "- %s.comparator_coverage_heatmap.pdf\n" "${figure_prefix}"
        printf -- "- %s.candidate_ont_support.pdf\n" "${figure_prefix}"
    } > "${hap_dir}/review_checklist.txt"
done

printf "haplotype\tcandidate_summary\tdepth_summary\tont_support_summary\tont_supported_candidates\tself_redundancy_summary\n" > "${OUT_DIR}/manifest.tsv"
for hap in "${HAPS[@]}"; do
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${hap}" \
        "${OUT_DIR}/${hap}/comparative_alignment/pzi_${hap}.non_aligned_candidate_summary.tsv" \
        "${OUT_DIR}/${hap}/ont_depth/pzi_${hap}.ont.primary.mapq${MIN_MAPQ}.depth_by_chromosome.tsv" \
        "${OUT_DIR}/${hap}/candidate_support/pzi_${hap}.candidate_ont_support_summary.tsv" \
        "${OUT_DIR}/${hap}/candidate_support/pzi_${hap}.candidate_ont_supported.tsv" \
        "${OUT_DIR}/${hap}/self_redundancy/pzi_${hap}.cross_chromosome_alignment_summary.tsv" \
        >> "${OUT_DIR}/manifest.tsv"
done

echo
echo "=== PZI expansion ONT support QC completed ==="
echo "Manifest: ${OUT_DIR}/manifest.tsv"
