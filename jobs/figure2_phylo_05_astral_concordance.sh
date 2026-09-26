#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
ALIGN="${BASE}/04_alignments"
GENES="${BASE}/05_gene_trees"
OUT="${BASE}/06_species_trees"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-8}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ASTRAL_JAR="${ASTRAL_JAR:-${PROJECT_ROOT}/ASTRAL-5.7.1/Astral/astral.5.7.1.jar}"
ASTRAL_MEM="${ASTRAL_MEM:-64G}"
ASTRAL_BOOTSTRAPS="${ASTRAL_BOOTSTRAPS:-1000}"
SEED="${SEED:-20260813}"

test -e "${GENES}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_04 first." >&2; exit 3; }
test -s "${ASTRAL_JAR}" || { echo "[ERROR] Missing ASTRAL jar: ${ASTRAL_JAR}" >&2; exit 2; }
test ! -e "${OUT}" || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 3; }
mkdir -p "${OUT}"

set +u
source "${CONDA_SH}"
conda activate "${TREE_ENV}"
set -u
IQTREE_BIN=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then IQTREE_BIN="$(command -v "${candidate}")"; break; fi
done
test -n "${IQTREE_BIN}" || { echo "[ERROR] IQ-TREE unavailable" >&2; exit 2; }
command -v java >/dev/null || { echo "[ERROR] java unavailable" >&2; exit 2; }

for dataset in main sensitivity; do
    DOUT="${OUT}/${dataset}"
    mkdir -p "${DOUT}/astral"
    find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | \
        while read -r name; do cat "${GENES}/${dataset}/genes/${name}"; done > "${DOUT}/all_gene_trees.uncontracted.tre"

    for threshold in 10 30 50; do
        "${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_contract_trees.py" \
            --tree-dir "${GENES}/${dataset}/genes" \
            --threshold "${threshold}" \
            --out "${DOUT}/all_gene_trees.bs${threshold}.tre"
        java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
            -i "${DOUT}/all_gene_trees.bs${threshold}.tre" \
            -t 2 \
            -o "${DOUT}/astral/species_tree.bs${threshold}.full_annotation.tre" \
            2> "${DOUT}/astral/species_tree.bs${threshold}.log"
    done
    java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
        -i "${DOUT}/all_gene_trees.bs10.tre" \
        -t 0 \
        -o "${DOUT}/astral/species_tree.bs10.clean.tre" \
        2> "${DOUT}/astral/species_tree.bs10.clean.log"

    find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | \
        sed 's/\.treefile$/.ufboot/' | \
        awk -v dir="${GENES}/${dataset}/genes" '{print dir "/" $0}' > "${DOUT}/astral/bootstrap_file_list.txt"

    java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
        -i "${DOUT}/all_gene_trees.bs10.tre" \
        -b "${DOUT}/astral/bootstrap_file_list.txt" \
        -g -r "${ASTRAL_BOOTSTRAPS}" -s "${SEED}" \
        -o "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" \
        2> "${DOUT}/astral/gene_site_bootstrap.log"
    tail -n 1 "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" > "${DOUT}/astral/gene_site_bootstrap.final_tree.tre"

    # ASTRAL 5.7.1 emits q1/q2/q3, quartet frequencies and posterior
    # annotations in the -t 2 full-annotation tree. It does not guarantee a
    # separate freqQuad.csv file, so retain a clearly named q123 copy instead.
    test -s "${DOUT}/astral/species_tree.bs10.full_annotation.tre" || {
        echo "[ERROR] Missing ASTRAL full-annotation tree" >&2
        exit 3
    }
    cp "${DOUT}/astral/species_tree.bs10.full_annotation.tre" \
       "${DOUT}/astral/species_tree.bs10.q123.tre"
    printf '%s\n' \
        "q1/q2/q3 and quartet-frequency annotations are stored in species_tree.bs10.full_annotation.tre (ASTRAL -t 2)." \
        > "${DOUT}/astral/species_tree.bs10.q123.log"

    "${IQTREE_BIN}" \
        -s "${ALIGN}/${dataset}/qc/concat_protein.fa" \
        -p "${ALIGN}/${dataset}/qc/concat_partitions.nex" \
        -m MFP+MERGE \
        -B 1000 --alrt 1000 --bnni \
        -T "${THREADS}" -seed "${SEED}" \
        --prefix "${DOUT}/concat_partitioned"

    "${IQTREE_BIN}" \
        -t "${DOUT}/astral/species_tree.bs10.clean.tre" \
        --gcf "${DOUT}/all_gene_trees.bs10.tre" \
        -s "${ALIGN}/${dataset}/qc/concat_protein.fa" \
        --scfl 1000 \
        -T "${THREADS}" -seed "${SEED}" \
        --prefix "${DOUT}/concordance_lscf"
done

{
    "${IQTREE_BIN}" --version
    java -version
    echo "ASTRAL_JAR=${ASTRAL_JAR}"
    sha256sum "${ASTRAL_JAR}"
} > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
