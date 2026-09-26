#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PACKAGE_DIR="${HAP_RENAME_DIR}"
GENOME_DIR="${PACKAGE_DIR}/genome"
GFF_DIR="${PACKAGE_DIR}/gff"
CDS_DIR="${PACKAGE_DIR}/cds"
PEP_DIR="${PACKAGE_DIR}/pep"
EVIANN_ROOT="${EVIANN_ROOT:-${ROOT_DIR}/annotation/genes/final_haplotypes_eviann}"
WORK_DIR="${WORK_DIR:-${ROOT_DIR}/annotation/genes/primary_isoform_busco}"
PRIMARY_DIR="${WORK_DIR}/primary"
BUSCO_DIR="${WORK_DIR}/busco"
SUMMARY_DIR="${WORK_DIR}/summary"
THREADS="${THREADS:-24}"
FORCE="${FORCE:-0}"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
SAVED_PACKAGE_CHECKSUMS="${WORK_DIR}/.haploid_rename.checksums.md5"
FORCE_EFFECTIVE="${FORCE}"
RUN_BUSCO="${RUN_BUSCO:-1}"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
PY_ENV="${PY_ENV:-haphic}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
SCRIPT="${ROOT_DIR}/scripts/eviann_primary_isoform_qc.py"

LABELS=(hapA1 hapA2 hapB1 hapB2 A1A2 B1B2 A1B1 A1B2 A2B1 A2B2 A1A2B1B2)

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "primary isoform Python script"
require_file "${BUSCO_LINEAGE}" "BUSCO lineage dataset"
mkdir -p "${PRIMARY_DIR}" "${BUSCO_DIR}" "${SUMMARY_DIR}" "${ROOT_DIR}/logs"
require_file "${PACKAGE_CHECKSUMS}" "renamed haplotype package checksums"

if ! cmp -s "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"; then
    FORCE_EFFECTIVE=1
    echo "[INFO] Renamed haplotype package changed; rebuilding primary sets and BUSCO results."
fi

for hap in hapA1 hapA2 hapB1 hapB2; do
    require_file "${GENOME_DIR}/pzi_${hap}.fa" "${hap} genome FASTA"
    require_file "${GFF_DIR}/pzi_${hap}.gff3" "${hap} EviAnn GFF3"
    require_file "${CDS_DIR}/pzi_${hap}.cds.fa" "${hap} EviAnn CDS"
    require_file "${PEP_DIR}/pzi_${hap}.pep.fa" "${hap} EviAnn protein"
    require_file "${EVIANN_ROOT}/${hap}/postprocess/pzi_${hap}.eviann.gene_model_qc.tsv" "${hap} EviAnn gene model QC"
done

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Build primary isoform sets from EviAnn annotation ==="
echo "Package : ${HAP_RENAME_DIR}"
echo "Genome  : ${GENOME_DIR}"
echo "GFF     : ${GFF_DIR}"
echo "CDS     : ${CDS_DIR}"
echo "Protein : ${PEP_DIR}"
echo "EviAnn  : ${EVIANN_ROOT}"
echo "Primary : ${PRIMARY_DIR}"
echo "BUSCO   : ${RUN_BUSCO}"
echo "Labels  : ${LABELS[*]}"
echo

if [ ! -s "${PRIMARY_DIR}/primary_summary.tsv" ] || [ "${FORCE_EFFECTIVE}" = "1" ]; then
    "${PYTHON_BIN}" "${SCRIPT}" build-primary \
        --package-dir "${PACKAGE_DIR}" \
        --eviann-root "${EVIANN_ROOT}" \
        --out-dir "${PRIMARY_DIR}"
else
    echo "Existing primary summary found, skipping primary rebuild: ${PRIMARY_DIR}/primary_summary.tsv"
fi

for label in "${LABELS[@]}"; do
    if [[ "${label}" == hap* ]]; then
        pep="${PRIMARY_DIR}/pzi_${label}.primary.pep.fa"
    else
        pep="${PRIMARY_DIR}/pzi_${label}.primary.pep.fa"
    fi
    require_file "${pep}" "${label} primary protein FASTA"
done

if [ "${RUN_BUSCO}" = "1" ]; then
    echo "=== BUSCO protein-mode for primary isoform sets ==="
    conda activate "${BUSCO_ENV}"
    for label in "${LABELS[@]}"; do
        pep="${PRIMARY_DIR}/pzi_${label}.primary.pep.fa"
        run_name="pzi_${label}_primary_protein_eudicots"
        existing_summary=$(find "${BUSCO_DIR}/${run_name}" -name "short_summary*.txt" -type f 2>/dev/null | head -n 1 || true)
        if [ -z "${existing_summary}" ] || [ "${FORCE_EFFECTIVE}" = "1" ]; then
            cd "${BUSCO_DIR}"
            busco \
                -i "${pep}" \
                -m proteins \
                -l "${BUSCO_LINEAGE}" \
                -c "${THREADS}" \
                --offline \
                -o "${run_name}" \
                --force
        else
            echo "Existing BUSCO result found for ${label}: ${existing_summary}"
        fi
    done
fi

conda activate "${PY_ENV}"
"${PYTHON_BIN}" "${SCRIPT}" summarize-busco \
    --busco-dir "${BUSCO_DIR}" \
    --out "${SUMMARY_DIR}/primary_busco_summary.tsv"

cp -f "${PRIMARY_DIR}/primary_summary.tsv" "${SUMMARY_DIR}/primary_summary.tsv"
cp -f "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"

cat > "${SUMMARY_DIR}/README.md" <<'EOF'
Primary isoform BUSCO comparison.

Inputs include four single haplotypes, same-group combinations, all A+B combinations,
and the full A1+A2+B1+B2 set. BUSCO was run in protein mode against eudicotyledons_odb12.

Interpretation:
- Complete = Single-copy + Duplicated BUSCO percentage.
- Single-copy = BUSCO genes found once in the input protein set.
- Duplicated = BUSCO genes found multiple times; expected to increase in A+B and all-haplotype combinations.
- Fragmented = partial BUSCO matches.
- Missing = BUSCO genes not detected.

Primary isoform inputs are better than full isoform inputs for BUSCO, OrthoFinder, RBH,
and gene-family analyses because each gene contributes only one representative protein.
EOF

echo "=== Primary isoform BUSCO workflow finished ==="
echo "Primary files : ${PRIMARY_DIR}"
echo "Summary files : ${SUMMARY_DIR}"
echo "BUSCO summary : ${SUMMARY_DIR}/primary_busco_summary.tsv"
