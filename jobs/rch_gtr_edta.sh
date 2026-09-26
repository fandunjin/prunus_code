#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
THREADS="${THREADS:-24}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
TARGET="${TARGET:-both}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

case "${TARGET}" in
    gtr|rch|both) ;;
    *) die "Unsupported TARGET=${TARGET} (use gtr, rch or both)" ;;
esac

require_file "${CONDA_SH}" "conda initialization"
if [[ ! -s "${EDTA_CDS}" && -s "${EDTA_CDS}.gz" ]]; then
    EDTA_CDS="${EDTA_CDS}.gz"
fi
require_file "${EDTA_CDS}" "EDTA CDS filter"

run_target() {
    local target="$1"
    local input out_dir work_fasta

    case "${target}" in
        gtr)
            input="${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna"
            out_dir="${TE_ROOT}/gtr"
            work_fasta="${out_dir}/gtr.fa"
            ;;
        rch)
            input="${ROOT_DIR}/data/rch/GCF_002994745.2_RchiOBHm-V2_genomic.fna.gz"
            out_dir="${TE_ROOT}/rch"
            work_fasta="${out_dir}/rch.fa"
            ;;
    esac

    require_file "${input}" "${target} genome input"
    mkdir -p "${out_dir}" "${ROOT_DIR}/logs"

    if [[ ! -s "${work_fasta}" ]]; then
        case "${input}" in
            *.gz) gzip -dc "${input}" > "${work_fasta}" ;;
            *)    cp -p "${input}" "${work_fasta}" ;;
        esac
    fi
    require_file "${work_fasta}" "${target} EDTA work FASTA"

    cd "${out_dir}"
    echo "=== Independent EDTA: ${target} ==="
    echo "Input : ${input}"
    echo "Work  : ${work_fasta}"
    echo "Out   : ${out_dir}"
    echo "CDS   : ${EDTA_CDS}"
    echo "Mode  : EDTA.pl --step all --species others --anno 1 --sensitive 1"

    if [[ ! -s "${work_fasta}.mod.EDTA.TEanno.gff3" ]]; then
        cmd=(EDTA.pl --genome "${work_fasta}" --species others --step all --anno 1 --sensitive 1 --force "${EDTA_FORCE}" --overwrite "${EDTA_OVERWRITE}" --threads "${THREADS}")
        if [[ -s "${EDTA_CDS}" ]]; then
            cmd+=(--cds "${EDTA_CDS}")
        fi
        "${cmd[@]}"
    else
        echo "[SKIP] Existing EDTA result: ${work_fasta}.mod.EDTA.TEanno.gff3"
    fi

    for output in \
        "${work_fasta}.mod.EDTA.TEanno.gff3" \
        "${work_fasta}.mod.EDTA.TEanno.sum" \
        "${work_fasta}.mod.EDTA.TElib.fa" \
        "${work_fasta}.mod.MAKER.masked"; do
        require_file "${output}" "${target} EDTA output"
    done

    printf 'target\tgenome\tte_gff3\tte_summary\tte_library\tmasked_genome\tstatus\n' > "${out_dir}/edta.completed.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\tpass\n' "${target}" "${work_fasta}" "${work_fasta}.mod.EDTA.TEanno.gff3" "${work_fasta}.mod.EDTA.TEanno.sum" "${work_fasta}.mod.EDTA.TElib.fa" "${work_fasta}.mod.MAKER.masked" >> "${out_dir}/edta.completed.tsv"
    echo "[OK] ${target} EDTA outputs: ${out_dir}"
}

source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
command -v EDTA.pl >/dev/null 2>&1 || die "EDTA.pl not found in environment ${EDTA_ENV}"

case "${TARGET}" in
    gtr|rch) run_target "${TARGET}" ;;
    both)    run_target gtr; run_target rch ;;
esac

echo "[DONE] EDTA outgroup runs finished."
