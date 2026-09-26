#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PYTHON_ENV="${PYTHON_ENV:-orthofinder}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/extract_swissprot_metadata.py}"
FUNCTIONAL_ROOT="${FUNCTIONAL_ROOT:-${ROOT_DIR}/annotation/functional}"
SWISSPROT_DIR="${SWISSPROT_DIR:-${HOME}/CX/database/SwissProt}"
SWISSPROT_DAT="${SWISSPROT_DAT:-${SWISSPROT_DIR}/uniprot_sprot.dat.gz}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "Swiss-Prot metadata parser"
require_file "${SWISSPROT_DAT}" "Swiss-Prot flat-file metadata"
[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || die "FORCE must be 0 or 1"

source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "${PYTHON_BIN} not found in ${PYTHON_ENV}"
mkdir -p "${FUNCTIONAL_ROOT}" "${ROOT_DIR}/logs"

protein_result() {
    case "$1" in
        pzi) printf '%s' "${FUNCTIONAL_ROOT}/pzi/swissprot.diamond.tsv" ;;
        psp) printf '%s' "${FUNCTIONAL_ROOT}/psp/swissprot.diamond.tsv" ;;
        ppa) printf '%s' "${FUNCTIONAL_ROOT}/ppa/swissprot.diamond.tsv" ;;
        *) die "Unknown species: $1 (choose pzi, psp or ppa)" ;;
    esac
}

run_species() {
    local species="$1"
    local diamond
    local out_dir="${FUNCTIONAL_ROOT}/${species}"
    diamond="$(protein_result "${species}")"
    require_file "${diamond}" "${species} Swiss-Prot DIAMOND result"
    mkdir -p "${out_dir}"

    if [[ "${FORCE}" == 0 && -s "${out_dir}/swissprot_go.tsv" && \
          -s "${out_dir}/swissprot_kegg.tsv" && \
          -s "${out_dir}/swissprot_ko.tsv" && \
          -s "${out_dir}/swissprot_ec.tsv" && \
          -s "${out_dir}/swissprot_interpro.tsv" && \
          -s "${out_dir}/swissprot_functional_annotation.tsv" ]]; then
        echo "[SKIP] ${species}: all Swiss-Prot metadata outputs already exist"
        return 0
    fi

    "${PYTHON_BIN}" "${SCRIPT}" \
        --species "${species}" \
        --diamond "${diamond}" \
        --swissprot-dat "${SWISSPROT_DAT}" \
        --out-dir "${out_dir}"
    for file in swissprot_go.tsv swissprot_kegg.tsv swissprot_ko.tsv \
                swissprot_ec.tsv swissprot_interpro.tsv \
                swissprot_functional_annotation.tsv; do
        require_file "${out_dir}/${file}" "${species} ${file}"
    done
    echo "[OK] ${species}: Swiss-Prot metadata transfer complete"
}

read -r -a SPECIES_LIST <<< "${SPECIES_SET//,/ }"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"
for species in "${SPECIES_LIST[@]}"; do
    run_species "${species}"
done

{
    echo "Swiss-Prot metadata transfer manifest"
    echo "Species requested: ${SPECIES_SET}"
    echo "Swiss-Prot metadata: ${SWISSPROT_DAT}"
    echo "DIAMOND results: ${FUNCTIONAL_ROOT}/<species>/swissprot.diamond.tsv"
    echo "Outputs: swissprot_go.tsv, swissprot_kegg.tsv, swissprot_ko.tsv, swissprot_ec.tsv, swissprot_interpro.tsv, swissprot_functional_annotation.tsv"
    echo "Method: accession join to the same Swiss-Prot release; homology transfer, not de novo prediction"
    echo "GO source: Swiss-Prot DR GO records"
    echo "KEGG/KO/EC/InterPro source: Swiss-Prot cross-reference records"
} > "${FUNCTIONAL_ROOT}/swissprot_metadata_transfer_manifest.txt"
echo "[OK] Swiss-Prot metadata extraction completed under ${FUNCTIONAL_ROOT}"
