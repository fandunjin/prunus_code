#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PACKAGE_DIR="${HAP_RENAME_DIR}"
GENOME_DIR="${PACKAGE_DIR}/genome"
GFF_DIR="${PACKAGE_DIR}/gff"
CDS_DIR="${PACKAGE_DIR}/cds"
PEP_DIR="${PACKAGE_DIR}/pep"
WORK_DIR="${WORK_DIR:-${ROOT_DIR}/haploid/internal_compare_eviann}"
OUT_DIR="${WORK_DIR}/summary"
SYN_DIR="${WORK_DIR}/synteny"
RBH_DIR="${WORK_DIR}/protein_rbh"
BUSCO_DIR="${WORK_DIR}/busco"
THREADS="${THREADS:-24}"
MINIMAP_PRESET="${MINIMAP_PRESET:-asm10}"
RUN_BUSCO="${RUN_BUSCO:-1}"
FORCE="${FORCE:-0}"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
SAVED_PACKAGE_CHECKSUMS="${WORK_DIR}/.haploid_rename.checksums.md5"
FORCE_EFFECTIVE="${FORCE}"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
PY_ENV="${PY_ENV:-haphic}"
RBH_ENV="${RBH_ENV:-eviann}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
SCRIPT="${ROOT_DIR}/scripts/subgenome_internal_compare.py"
DOTPLOT_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py"
HEATMAP_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_heatmap.py"

HAPS=(hapA1 hapA2 hapB1 hapB2)
SYNTENY_PAIRS=(
    "hapA1 hapA2 A1_vs_A2"
    "hapB1 hapB2 B1_vs_B2"
    "hapA1 hapB1 A1_vs_B1"
    "hapA1 hapB2 A1_vs_B2"
    "hapA2 hapB1 A2_vs_B1"
    "hapA2 hapB2 A2_vs_B2"
)
RBH_PAIRS=(
    "hapA1 hapA2"
    "hapB1 hapB2"
    "hapA1 hapB1"
    "hapA1 hapB2"
    "hapA2 hapB1"
    "hapA2 hapB2"
)

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

remove_if_forced() {
    if [ "${FORCE_EFFECTIVE}" = "1" ]; then
        rm -f "$@"
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "subgenome internal comparison Python script"
require_file "${DOTPLOT_SCRIPT}" "minimap dotplot script"
require_file "${HEATMAP_SCRIPT}" "minimap heatmap script"
mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${SYN_DIR}" "${RBH_DIR}" "${BUSCO_DIR}" "${ROOT_DIR}/logs"
require_file "${PACKAGE_CHECKSUMS}" "renamed haplotype package checksums"

if ! cmp -s "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"; then
    FORCE_EFFECTIVE=1
    echo "[INFO] Renamed haplotype package changed; rebuilding cached comparisons."
fi

for hap in "${HAPS[@]}"; do
    require_file "${GENOME_DIR}/pzi_${hap}.fa" "${hap} genome FASTA"
    require_file "${GFF_DIR}/pzi_${hap}.gff3" "${hap} GFF3"
    require_file "${CDS_DIR}/pzi_${hap}.cds.fa" "${hap} CDS FASTA"
    require_file "${PEP_DIR}/pzi_${hap}.pep.fa" "${hap} protein FASTA"
done

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== PZI final haplotype/subgenome internal comparison ==="
echo "Package : ${HAP_RENAME_DIR}"
echo "Genome  : ${GENOME_DIR}"
echo "GFF     : ${GFF_DIR}"
echo "CDS     : ${CDS_DIR}"
echo "Protein : ${PEP_DIR}"
echo "Work    : ${WORK_DIR}"
echo "Out     : ${OUT_DIR}"
echo "Threads : ${THREADS}"
echo "Preset  : ${MINIMAP_PRESET}"
echo "BUSCO   : ${RUN_BUSCO}"
echo "RBH env : ${RBH_ENV}"
echo "Note    : using chromosome-renamed four-haplotype package."
echo

echo ">>> Genome synteny: A1 vs A2, B1 vs B2, and A-vs-B pairs"
for item in "${SYNTENY_PAIRS[@]}"; do
    set -- ${item}
    ref_hap="$1"
    qry_hap="$2"
    label="$3"
    ref_fa="${GENOME_DIR}/pzi_${ref_hap}.fa"
    qry_fa="${GENOME_DIR}/pzi_${qry_hap}.fa"
    paf="${SYN_DIR}/${label}.${MINIMAP_PRESET}.paf"
    summary="${SYN_DIR}/${label}.summary.tsv"
    dotplot="${SYN_DIR}/${label}.${MINIMAP_PRESET}.dotplot.svg"
    heatmap="${SYN_DIR}/${label}.${MINIMAP_PRESET}.heatmap.svg"
    remove_if_forced "${paf}" "${summary}" "${dotplot}" "${heatmap}"
    if [ ! -s "${paf}" ]; then
        minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${ref_fa}" "${qry_fa}" > "${paf}"
    fi
    if [ ! -s "${summary}" ]; then
        "${PYTHON_BIN}" "${SCRIPT}" paf-summary --paf "${paf}" --out "${summary}"
    fi
    if [ ! -s "${dotplot}" ]; then
        "${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
            --paf "${paf}" \
            --out-svg "${dotplot}" \
            --title "${label} ${MINIMAP_PRESET} synteny" \
            --min-alignment-len 50000 \
            --min-query-aggregate 100000 \
            --min-target-aggregate 100000 \
            --min-visible-px 1.2
    fi
    if [ ! -s "${heatmap}" ]; then
        "${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
            --summary "${summary}" \
            --out-svg "${heatmap}" \
            --title "${label} chromosome-pair aligned bp" \
            --x-label "${qry_hap}" \
            --y-label "${ref_hap}"
    fi
done

echo ">>> Protein reciprocal best-hit analysis"
if ! command -v diamond >/dev/null 2>&1 && ! { command -v makeblastdb >/dev/null 2>&1 && command -v blastp >/dev/null 2>&1; }; then
    echo "Protein alignment tools were not found in ${PY_ENV}; trying RBH_ENV=${RBH_ENV}."
    conda activate "${RBH_ENV}"
fi
if command -v diamond >/dev/null 2>&1; then
    for hap in "${HAPS[@]}"; do
        db="${RBH_DIR}/pzi_${hap}.dmnd"
        if [ ! -s "${db}" ] || [ "${FORCE_EFFECTIVE}" = "1" ]; then
            if [ "${FORCE_EFFECTIVE}" = "1" ]; then
                rm -f "${db}"
            fi
            diamond makedb --in "${PEP_DIR}/pzi_${hap}.pep.fa" --db "${RBH_DIR}/pzi_${hap}" --threads "${THREADS}"
        fi
    done
    for item in "${RBH_PAIRS[@]}"; do
        set -- ${item}
        q="$1"
        t="$2"
        out1="${RBH_DIR}/${q}_vs_${t}.blastp.tsv"
        out2="${RBH_DIR}/${t}_vs_${q}.blastp.tsv"
        remove_if_forced "${out1}" "${out2}"
        if [ ! -s "${out1}" ]; then
            diamond blastp \
                --query "${PEP_DIR}/pzi_${q}.pep.fa" \
                --db "${RBH_DIR}/pzi_${t}.dmnd" \
                --out "${out1}" \
                --outfmt 6 qseqid sseqid pident length evalue bitscore \
                --max-target-seqs 5 \
                --evalue 1e-5 \
                --threads "${THREADS}"
        fi
        if [ ! -s "${out2}" ]; then
            diamond blastp \
                --query "${PEP_DIR}/pzi_${t}.pep.fa" \
                --db "${RBH_DIR}/pzi_${q}.dmnd" \
                --out "${out2}" \
                --outfmt 6 qseqid sseqid pident length evalue bitscore \
                --max-target-seqs 5 \
                --evalue 1e-5 \
                --threads "${THREADS}"
        fi
    done
elif command -v makeblastdb >/dev/null 2>&1 && command -v blastp >/dev/null 2>&1; then
    for hap in "${HAPS[@]}"; do
        db_prefix="${RBH_DIR}/pzi_${hap}.blastdb"
        if [ ! -s "${db_prefix}.pin" ] || [ "${FORCE_EFFECTIVE}" = "1" ]; then
            if [ "${FORCE_EFFECTIVE}" = "1" ]; then
                rm -f "${db_prefix}".*
            fi
            makeblastdb -in "${PEP_DIR}/pzi_${hap}.pep.fa" -dbtype prot -out "${db_prefix}"
        fi
    done
    for item in "${RBH_PAIRS[@]}"; do
        set -- ${item}
        q="$1"
        t="$2"
        out1="${RBH_DIR}/${q}_vs_${t}.blastp.tsv"
        out2="${RBH_DIR}/${t}_vs_${q}.blastp.tsv"
        remove_if_forced "${out1}" "${out2}"
        if [ ! -s "${out1}" ]; then
            blastp \
                -query "${PEP_DIR}/pzi_${q}.pep.fa" \
                -db "${RBH_DIR}/pzi_${t}.blastdb" \
                -out "${out1}" \
                -outfmt "6 qseqid sseqid pident length evalue bitscore" \
                -max_target_seqs 5 \
                -evalue 1e-5 \
                -num_threads "${THREADS}"
        fi
        if [ ! -s "${out2}" ]; then
            blastp \
                -query "${PEP_DIR}/pzi_${t}.pep.fa" \
                -db "${RBH_DIR}/pzi_${q}.blastdb" \
                -out "${out2}" \
                -outfmt "6 qseqid sseqid pident length evalue bitscore" \
                -max_target_seqs 5 \
                -evalue 1e-5 \
                -num_threads "${THREADS}"
        fi
    done
else
    echo "ERROR: neither diamond nor BLAST+ was found in PATH." >&2
    echo "Tried PY_ENV=${PY_ENV} and RBH_ENV=${RBH_ENV}. Set RBH_ENV to an environment containing diamond or BLAST+." >&2
    exit 1
fi

if [ "${RUN_BUSCO}" = "1" ]; then
    echo ">>> BUSCO protein-mode comparison for four haplotypes"
    conda activate "${BUSCO_ENV}"
    require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
    for hap in "${HAPS[@]}"; do
        run_name="pzi_${hap}_protein_eudicots"
        existing_busco_summary=$(find "${BUSCO_DIR}/${run_name}" -name "short_summary*.txt" -type f 2>/dev/null | head -n 1 || true)
        if [ -z "${existing_busco_summary}" ] || [ "${FORCE_EFFECTIVE}" = "1" ]; then
            cd "${BUSCO_DIR}"
            busco \
                -i "${PEP_DIR}/pzi_${hap}.pep.fa" \
                -m proteins \
                -l "${BUSCO_LINEAGE}" \
                -c "${THREADS}" \
                --offline \
                -o "${run_name}" \
                --force
        fi
    done
    conda activate "${PY_ENV}"
fi

echo ">>> Collect comparison summaries and recommended downstream representative"
conda activate "${PY_ENV}"
"${PYTHON_BIN}" "${SCRIPT}" collect \
    --package-dir "${PACKAGE_DIR}" \
    --work-dir "${WORK_DIR}" \
    --out-dir "${OUT_DIR}"

cp -f "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"

echo "=== Internal comparison completed ==="
echo "Summary dir:"
echo "${OUT_DIR}"
echo
echo "Key outputs:"
echo "${OUT_DIR}/synteny_summary.tsv"
echo "${OUT_DIR}/hap_gene_count_summary.tsv"
echo "${OUT_DIR}/AB_gene_number_difference.tsv"
echo "${OUT_DIR}/rbh_gene_retention_summary.tsv"
echo "${OUT_DIR}/gene_loss_candidates.tsv"
echo "${OUT_DIR}/busco_summary.tsv"
echo "${OUT_DIR}/hap_representative_ranking.tsv"
echo "${OUT_DIR}/downstream_representative_recommendation.tsv"
echo "${OUT_DIR}/downstream_representative_AplusB/"
