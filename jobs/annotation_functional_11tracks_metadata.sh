#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# No DIAMOND search or gene prediction is rerun.
# The --swissprot-dat flat file must be the SAME Swiss-Prot release used to
# build annotation/functional/work/swissprot.dmnd.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-orthofinder}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/extract_swissprot_metadata.py}"
FUNCTIONAL_ROOT="${FUNCTIONAL_ROOT:-${ROOT_DIR}/annotation/functional}"
SWISSPROT_DIR="${SWISSPROT_DIR:-${HOME}/CX/database/SwissProt}"
SWISSPROT_DAT="${SWISSPROT_DAT:-${SWISSPROT_DIR}/uniprot_sprot.dat.gz}"
SPECIES_SET="${SPECIES_SET:-pcam pav pspe pmal par pmu pdu psal ppe gtr rch}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || die "FORCE must be 0 or 1"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "Swiss-Prot metadata parser"
require_file "${SWISSPROT_DAT}" "Swiss-Prot flat-file metadata"

source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "${PYTHON_BIN} not found in ${PYTHON_ENV}"

mkdir -p "${FUNCTIONAL_ROOT}" "${ROOT_DIR}/logs"

run_track() {
    local track="$1"
    local out_dir="${FUNCTIONAL_ROOT}/${track}"
    local diamond="${out_dir}/swissprot.diamond.tsv"

    case "$track" in
        pzi|psp|ppa)
            die "refusing to annotate ${track}: racemose metadata is owned by jobs/annotation_functional_swissprot_metadata.sh" ;;
    esac

    require_file "${diamond}" "${track} Swiss-Prot DIAMOND result (run jobs/annotation_functional_11tracks.sh first)"
    mkdir -p "${out_dir}"

    if [[ "${FORCE}" == 0 && -s "${out_dir}/swissprot_go.tsv" && \
          -s "${out_dir}/swissprot_kegg.tsv" && \
          -s "${out_dir}/swissprot_ko.tsv" && \
          -s "${out_dir}/swissprot_ec.tsv" && \
          -s "${out_dir}/swissprot_interpro.tsv" && \
          -s "${out_dir}/swissprot_functional_annotation.tsv" ]]; then
        echo "[SKIP] ${track}: all Swiss-Prot metadata outputs already exist"
        return 0
    fi

    "${PYTHON_BIN}" "${SCRIPT}" \
        --species "${track}" \
        --diamond "${diamond}" \
        --swissprot-dat "${SWISSPROT_DAT}" \
        --out-dir "${out_dir}"

    for file in swissprot_go.tsv swissprot_kegg.tsv swissprot_ko.tsv \
                swissprot_ec.tsv swissprot_interpro.tsv \
                swissprot_functional_annotation.tsv; do
        require_file "${out_dir}/${file}" "${track} ${file}"
    done
    echo "[OK] ${track}: Swiss-Prot metadata transfer complete"
}

read -r -a TRACKS <<< "${SPECIES_SET//,/ }"
[[ "${#TRACKS[@]}" -gt 0 ]] || die "SPECIES_SET is empty"

for track in "${TRACKS[@]}"; do
    run_track "${track}"
done

{
    echo "Swiss-Prot metadata transfer manifest (11 tracks)"
    echo "Tracks requested: ${SPECIES_SET}"
    echo "Swiss-Prot metadata: ${SWISSPROT_DAT}"
    echo "DIAMOND results: ${FUNCTIONAL_ROOT}/<track>/swissprot.diamond.tsv"
    echo "Outputs per track: swissprot_go.tsv, swissprot_kegg.tsv, swissprot_ko.tsv,"
    echo "                   swissprot_ec.tsv, swissprot_interpro.tsv,"
    echo "                   swissprot_functional_annotation.tsv"
    echo "Method: accession join to the same Swiss-Prot release; homology transfer"
    echo "GO source: Swiss-Prot DR GO records"
    echo "KEGG/KO/EC/InterPro source: Swiss-Prot cross-reference records"
} > "${FUNCTIONAL_ROOT}/swissprot_metadata_transfer_manifest_11tracks.txt"

echo "[OK] 11-track Swiss-Prot metadata extraction completed under ${FUNCTIONAL_ROOT}"
