#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
ALIGN="${BASE}/04_alignments"
GENES="${BASE}/05_gene_trees"
OUT="${BASE}/06_species_trees"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-48}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ASTRAL_JAR="${ASTRAL_JAR:-${PROJECT_ROOT}/ASTRAL-5.7.1/Astral/astral.5.7.1.jar}"
ASTRAL_MEM="${ASTRAL_MEM:-64G}"
ASTRAL_BOOTSTRAPS="${ASTRAL_BOOTSTRAPS:-1000}"
SEED="${SEED:-20260813}"

test -e "${GENES}/COMPLETE" || { echo "[ERROR] Missing 05_gene_trees/COMPLETE" >&2; exit 3; }
test -d "${OUT}" || { echo "[ERROR] Missing partial 06_species_trees directory" >&2; exit 3; }
test -s "${ASTRAL_JAR}" || { echo "[ERROR] Missing ASTRAL jar: ${ASTRAL_JAR}" >&2; exit 2; }

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

run_astral_dataset() {
    local dataset="$1"
    local DOUT="${OUT}/${dataset}"
    mkdir -p "${DOUT}/astral"

    find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | \
        while read -r name; do cat "${GENES}/${dataset}/genes/${name}"; done > "${DOUT}/all_gene_trees.uncontracted.tre"

    for threshold in 10 30 50; do
        if [ ! -s "${DOUT}/all_gene_trees.bs${threshold}.tre" ]; then
            "${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_contract_trees.py" \
                --tree-dir "${GENES}/${dataset}/genes" \
                --threshold "${threshold}" \
                --out "${DOUT}/all_gene_trees.bs${threshold}.tre"
        fi
        if [ ! -s "${DOUT}/astral/species_tree.bs${threshold}.full_annotation.tre" ]; then
            java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
                -i "${DOUT}/all_gene_trees.bs${threshold}.tre" \
                -t 2 \
                -o "${DOUT}/astral/species_tree.bs${threshold}.full_annotation.tre" \
                2> "${DOUT}/astral/species_tree.bs${threshold}.log"
        fi
    done

    if [ ! -s "${DOUT}/astral/species_tree.bs10.clean.tre" ]; then
        java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
            -i "${DOUT}/all_gene_trees.bs10.tre" \
            -t 0 \
            -o "${DOUT}/astral/species_tree.bs10.clean.tre" \
            2> "${DOUT}/astral/species_tree.bs10.clean.log"
    fi

    if [ ! -s "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" ]; then
        find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | \
            sed 's/\.treefile$/.ufboot/' | \
            awk -v dir="${GENES}/${dataset}/genes" '{print dir "/" $0}' > "${DOUT}/astral/bootstrap_file_list.txt"
        java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" \
            -i "${DOUT}/all_gene_trees.bs10.tre" \
            -b "${DOUT}/astral/bootstrap_file_list.txt" \
            -g -r "${ASTRAL_BOOTSTRAPS}" -s "${SEED}" \
            -o "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" \
            2> "${DOUT}/astral/gene_site_bootstrap.log"
    fi
    if [ ! -s "${DOUT}/astral/gene_site_bootstrap.final_tree.tre" ]; then
        tail -n 1 "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" > "${DOUT}/astral/gene_site_bootstrap.final_tree.tre"
    fi

    cp "${DOUT}/astral/species_tree.bs10.full_annotation.tre" \
       "${DOUT}/astral/species_tree.bs10.q123.tre"
    printf '%s\n' \
        "q1/q2/q3 and quartet-frequency annotations are stored in species_tree.bs10.full_annotation.tre (ASTRAL -t 2)." \
        > "${DOUT}/astral/species_tree.bs10.q123.log"

    if [ ! -s "${DOUT}/concat_partitioned.treefile" ]; then
        "${IQTREE_BIN}" \
            -s "${ALIGN}/${dataset}/qc/concat_protein.fa" \
            -p "${ALIGN}/${dataset}/qc/concat_partitions.nex" \
            -m MFP+MERGE \
            -B 1000 --alrt 1000 --bnni \
            -T "${THREADS}" -seed "${SEED}" \
            --prefix "${DOUT}/concat_partitioned"
    fi

    if [ ! -s "${DOUT}/concordance_lscf.cf.tree" ]; then
        "${IQTREE_BIN}" \
            -t "${DOUT}/astral/species_tree.bs10.clean.tre" \
            --gcf "${DOUT}/all_gene_trees.bs10.tre" \
            -s "${ALIGN}/${dataset}/qc/concat_protein.fa" \
            --scfl 1000 \
            -T "${THREADS}" -seed "${SEED}" \
            --prefix "${DOUT}/concordance_lscf"
    fi
}

run_astral_dataset main
run_astral_dataset sensitivity

{
    "${IQTREE_BIN}" --version
    java -version
    echo "ASTRAL_JAR=${ASTRAL_JAR}"
    sha256sum "${ASTRAL_JAR}"
} > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
