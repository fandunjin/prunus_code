#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
TREE_DIR="${ANALYSIS_DIR}/trees/main"
OVERWRITE="${OVERWRITE:-no}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
TREE_ENV="${TREE_ENV:-orthofinder}"
CODON_ENV="${CODON_ENV:-${TREE_ENV}}"
THREADS="${THREADS:-48}"
ASTRAL_JAR="${ASTRAL_JAR:-${PROJECT_ROOT}/ASTRAL-5.7.1/Astral/astral.5.7.1.jar}"
ASTRAL_MEM="${ASTRAL_MEM:-64G}"
IQTREE_BOOTSTRAPS="${IQTREE_BOOTSTRAPS:-1000}"
IQTREE_ALRT="${IQTREE_ALRT:-1000}"
RAXML_BOOTSTRAPS="${RAXML_BOOTSTRAPS:-1000}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ORTHOFINDER_LATEST="${ORTHOFINDER_LATEST:-${ANALYSIS_DIR}/orthofinder/results/main.latest.txt}"
SPECIES="${SPECIES:-pzi_A2B2,ppe,pmu_refseq,pdu,pcam,pspe,mdom,msyl,pcom,rch,fve}"

if [ -e "${TREE_DIR}" ] && [ "${OVERWRITE}" != "yes" ]; then
    echo "[ERROR] TREE_DIR already exists: ${TREE_DIR}" >&2
    echo "[ERROR] Move it aside or resubmit with OVERWRITE=yes after confirming it only contains failed/obsolete outputs." >&2
    exit 3
fi
if [ -e "${TREE_DIR}" ] && [ "${OVERWRITE}" = "yes" ]; then
    rm -rf "${TREE_DIR}"
fi
mkdir -p "${TREE_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${CODON_ENV}"
set -u

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN=python
    else
        echo "[ERROR] Neither ${PYTHON_BIN} nor python is available in CODON_ENV=${CODON_ENV}" >&2
        exit 2
    fi
fi

for exe in mafft trimal AMAS.py; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in CODON_ENV=${CODON_ENV}: ${exe}" >&2
        exit 2
    fi
done

if [[ ! -s "${ASTRAL_JAR}" ]]; then
    echo "[ERROR] ASTRAL_JAR is not set to an existing ASTRAL jar: ${ASTRAL_JAR}" >&2
    echo "[ERROR] Set ASTRAL_JAR=/path/to/astral.jar before submitting downstream_04." >&2
    exit 2
fi

if [ ! -s "${ORTHOFINDER_LATEST}" ]; then
    echo "[ERROR] Missing OrthoFinder latest pointer: ${ORTHOFINDER_LATEST}" >&2
    echo "[ERROR] Run downstream_03_orthofinder.sh successfully before downstream_04." >&2
    exit 3
fi
read -r of_dir < "${ORTHOFINDER_LATEST}"
if [ ! -d "${of_dir}" ]; then
    echo "[ERROR] OrthoFinder result directory from ${ORTHOFINDER_LATEST} does not exist: ${of_dir}" >&2
    exit 3
fi
orthogroups="${of_dir}/Orthogroups/Orthogroups.tsv"
if [ ! -s "${orthogroups}" ]; then
    echo "[ERROR] Missing OrthoFinder orthogroups table: ${orthogroups}" >&2
    exit 3
fi

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_extract_single_copy_for_pal2nal.py" \
    --canonical-manifest "${ANALYSIS_DIR}/canonical_manifest.tsv" \
    --orthogroups-tsv "${orthogroups}" \
    --species "${SPECIES}" \
    --out-dir "${TREE_DIR}"

if [ ! -s "${TREE_DIR}/single_copy_orthogroups.txt" ]; then
    echo "[ERROR] No usable single-copy orthogroups were extracted for codon-alignment construction." >&2
    echo "[ERROR] Check ${TREE_DIR}/single_copy_extract_stats.tsv for skipped orthogroups." >&2
    exit 3
fi

mkdir -p \
    "${TREE_DIR}/pal2nal_inputs/aa_alignments" \
    "${TREE_DIR}/pal2nal_inputs/aa_alignments_trimmed" \
    "${TREE_DIR}/codon_alignments" \
    "${TREE_DIR}/logs/mafft" \
    "${TREE_DIR}/logs/trimal"
for aa_fa in "${TREE_DIR}"/pal2nal_inputs/aa_unaligned/*.fa; do
    base="${aa_fa##*/}"
    base="${base%.fa}"
    cds_fa="${TREE_DIR}/pal2nal_inputs/cds/${base}.fa"
    aa_aln="${TREE_DIR}/pal2nal_inputs/aa_alignments/${base}.fa"
    trimmed_aa="${TREE_DIR}/pal2nal_inputs/aa_alignments_trimmed/${base}.fa"
    trimmed_codon="${TREE_DIR}/codon_alignments/${base}.fa"
    mafft_log="${TREE_DIR}/logs/mafft/${base}.log"
    trimal_aa_log="${TREE_DIR}/logs/trimal/${base}.aa.log"
    trimal_codon_log="${TREE_DIR}/logs/trimal/${base}.codon.log"

    if [ ! -s "${cds_fa}" ]; then
        echo "[ERROR] Missing CDS FASTA for codon-alignment input ${base}: ${cds_fa}" >&2
        exit 3
    fi

    mafft --auto --thread 1 "${aa_fa}" > "${aa_aln}" 2> "${mafft_log}"
    if [ ! -s "${aa_aln}" ]; then
        echo "[ERROR] MAFFT produced empty protein alignment for ${base}" >&2
        exit 3
    fi

    trimal -in "${aa_aln}" -out "${trimmed_aa}" -fasta -gt 0.8 -st 0.001 > "${trimal_aa_log}" 2>&1 || true
    if [ ! -s "${trimmed_aa}" ]; then
        echo "[WARN] trimAl produced empty protein alignment for ${base}; skipping this orthogroup." >&2
        rm -f "${trimmed_aa}" "${trimmed_codon}"
        echo -e "${base}\ttrimal_aa_empty" >> "${TREE_DIR}/codon_alignment_skipped.tsv"
        continue
    fi

    trimal -in "${aa_aln}" -backtrans "${cds_fa}" -out "${trimmed_codon}" -fasta -gt 0.8 -st 0.001 > "${trimal_codon_log}" 2>&1 || true
    if [ ! -s "${trimmed_codon}" ]; then
        echo "[WARN] trimAl -backtrans produced empty codon alignment for ${base}; skipping this orthogroup." >&2
        rm -f "${trimmed_codon}"
        echo -e "${base}\ttrimal_backtrans_empty" >> "${TREE_DIR}/codon_alignment_skipped.tsv"
        continue
    fi
done

codon_count=0
for codon_aln in "${TREE_DIR}"/codon_alignments/*.fa; do
    if [ -s "${codon_aln}" ]; then
        ((codon_count += 1))
    fi
done
if [ "${codon_count}" -lt 100 ]; then
    echo "[ERROR] Only ${codon_count} codon alignments were retained; this is too few for the publication-grade species tree workflow." >&2
    echo "[ERROR] Check ${TREE_DIR}/single_copy_extract_status_counts.tsv and ${TREE_DIR}/codon_alignment_skipped.tsv." >&2
    exit 3
fi

cd "${TREE_DIR}/codon_alignments"
AMAS.py concat \
    -i *.fa \
    -f fasta \
    -d dna \
    -u phylip \
    -t "${TREE_DIR}/concat_codon.phy" \
    -p "${TREE_DIR}/concat_partitions.txt" \
    -y raxml

if [ ! -s "${TREE_DIR}/concat_codon.phy" ] || [ ! -s "${TREE_DIR}/concat_partitions.txt" ]; then
    echo "[ERROR] Missing concatenated codon alignment or partition file in ${TREE_DIR}" >&2
    exit 3
fi

set +u
conda activate "${TREE_ENV}"
set -u

for exe in iqtree3 raxml-ng java; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in TREE_ENV=${TREE_ENV}: ${exe}" >&2
        exit 2
    fi
done

mkdir -p "${TREE_DIR}/gene_trees"
for aln in "${TREE_DIR}"/codon_alignments/*.fa; do
    base="${aln##*/}"
    base="${base%.fa}"
    iqtree3 \
        -s "${aln}" \
        -m MFP \
        -B "${IQTREE_BOOTSTRAPS}" \
        --alrt "${IQTREE_ALRT}" \
        -T 1 \
        --prefix "${TREE_DIR}/gene_trees/${base}"
done

cat "${TREE_DIR}"/gene_trees/*.treefile > "${TREE_DIR}/gene_trees/all_gene_trees.tre"
if [ ! -s "${TREE_DIR}/gene_trees/all_gene_trees.tre" ]; then
    echo "[ERROR] No gene tree files were generated under ${TREE_DIR}/gene_trees" >&2
    exit 3
fi

iqtree3 \
    -s "${TREE_DIR}/concat_codon.phy" \
    -p "${TREE_DIR}/concat_partitions.txt" \
    -m MFP \
    -B "${IQTREE_BOOTSTRAPS}" \
    --alrt "${IQTREE_ALRT}" \
    -T "${THREADS}" \
    --prefix "${TREE_DIR}/iqtree_concat"

raxml-ng \
    --all \
    --msa "${TREE_DIR}/concat_codon.phy" \
    --model "${TREE_DIR}/concat_partitions.txt" \
    --brlen scaled \
    --tree pars{10},rand{10} \
    --bs-trees "${RAXML_BOOTSTRAPS}" \
    --threads "${THREADS}" \
    --prefix "${TREE_DIR}/raxmlng_concat"

java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
    -i "${TREE_DIR}/gene_trees/all_gene_trees.tre" \
    -o "${TREE_DIR}/astral_species_tree.tre" \
    2> "${TREE_DIR}/astral_species_tree.log"

iqtree3 \
    -t "${TREE_DIR}/iqtree_concat.treefile" \
    --gcf "${TREE_DIR}/gene_trees/all_gene_trees.tre" \
    --prefix "${TREE_DIR}/concordance_gcf"

iqtree3 \
    -te "${TREE_DIR}/iqtree_concat.treefile" \
    -p "${TREE_DIR}/codon_alignments" \
    --scfl 100 \
    --prefix "${TREE_DIR}/concordance_scfl"
