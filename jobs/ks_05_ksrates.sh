#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
KAKS_ENV="${KAKS_ENV:-akr}"
JCVI_INPUT_ROOT="${JCVI_INPUT_ROOT:-${KS_ROOT}/jcvi/inputs}"
INPUT_ROOT="${INPUT_ROOT:-${KS_ROOT}/ksrates/input}"
RUN_ROOT="${RUN_ROOT:-${KS_ROOT}/ksrates/runs}"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_ksrates_inputs.py"
THREADS="${THREADS:-32}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
KSRATES_REVISION="${KSRATES_REVISION:-v2.1.1}"
KSRATES_PIPELINE="${KSRATES_PIPELINE:-VIB-PSB/ksrates}"
KSRATES_PROFILE="${KSRATES_PROFILE:-apptainer}"
KSRATES_ONLY="${KSRATES_ONLY:-pzi psp ppa}"
NEXTFLOW_CMD="${NEXTFLOW_CMD:-nextflow}"
APPTAINER_CMD="${APPTAINER_CMD:-apptainer}"
PREP_ONLY="${PREP_ONLY:-0}"
SKIP_PREP="${SKIP_PREP:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

mkdir -p "${KS_ROOT}/logs" "${INPUT_ROOT}" "${RUN_ROOT}"
rm -f "${KS_ROOT}/ksrates/COMPLETE"
require_file "${PREP_SCRIPT}"
[[ -e "${KS_ROOT}/jcvi/COMPLETE" ]] || die "Missing completion marker: ${KS_ROOT}/jcvi/COMPLETE"

set +u
source "${CONDA_SH}"
conda activate "${KAKS_ENV}"
set -u
command -v python >/dev/null 2>&1 || die "python unavailable"

if [[ "${SKIP_PREP}" == "0" ]]; then
    python "${PREP_SCRIPT}" \
        --jcvi-input-root "${JCVI_INPUT_ROOT}" \
        --outdir "${INPUT_ROOT}" \
        --bootstrap-iterations "${BOOTSTRAPS}" \
        --pipeline-version "${KSRATES_REVISION}" \
        --threads "${THREADS}"
else
    require_file "${INPUT_ROOT}/PREPARED"
fi

require_file "${INPUT_ROOT}/input_qc.tsv"
for focal in pzi psp ppa; do
    require_file "${INPUT_ROOT}/configs/config_${focal}.txt"
    require_file "${INPUT_ROOT}/sequences/${focal}.fasta"
    require_file "${INPUT_ROOT}/annotations/${focal}.gff3"
done
for reference in pcam ppe gtr; do
    require_file "${INPUT_ROOT}/sequences/${reference}.fasta"
done

if [[ "${PREP_ONLY}" == "1" ]]; then
    echo "[OK] PREP_ONLY=1; ksrates inputs are ready at ${INPUT_ROOT}"
    exit 0
fi

command -v "${NEXTFLOW_CMD}" >/dev/null 2>&1 || die "Nextflow unavailable"
if [[ "${KSRATES_PROFILE}" == "apptainer" ]]; then
    command -v "${APPTAINER_CMD}" >/dev/null 2>&1 || die "Apptainer unavailable"
elif [[ "${KSRATES_PROFILE}" == "singularity" ]]; then
    command -v singularity >/dev/null 2>&1 || die "Singularity unavailable"
else
    die "KSRATES_PROFILE must be apptainer or singularity"
fi

revision_args=()
if [[ ! -d "${KSRATES_PIPELINE}" ]]; then
    revision_args=(-r "${KSRATES_REVISION}")
fi

for focal in ${KSRATES_ONLY}; do
    case "${focal}" in
        pzi|psp|ppa) ;;
        *) die "Unsupported KSRATES_ONLY value: ${focal}" ;;
    esac
    run_dir="${RUN_ROOT}/${focal}"
    mkdir -p "${run_dir}"
    rm -f "${run_dir}/COMPLETE"
    cd "${run_dir}"
    "${NEXTFLOW_CMD}" run "${KSRATES_PIPELINE}" \
        "${revision_args[@]}" \
        -resume \
        -profile "${KSRATES_PROFILE}" \
        -c "${INPUT_ROOT}/configs/nextflow.config" \
        -work-dir "${run_dir}/work" \
        -with-report "${run_dir}/nextflow_report.html" \
        -with-timeline "${run_dir}/nextflow_timeline.html" \
        -with-trace "${run_dir}/nextflow_trace.tsv" \
        -with-dag "${run_dir}/nextflow_dag.html" \
        --config "${INPUT_ROOT}/configs/config_${focal}.txt" \
        --expert "${INPUT_ROOT}/configs/config_expert.txt"
    touch "${run_dir}/COMPLETE"
done

{
    echo "ksrates_revision=${KSRATES_REVISION}"
    echo "ksrates_pipeline=${KSRATES_PIPELINE}"
    echo "profile=${KSRATES_PROFILE}"
    echo "nextflow=$(${NEXTFLOW_CMD} -version 2>&1 | head -n 2 | tail -n 1)"
    echo "threads=${THREADS}"
    echo "ortholog_peak_bootstraps=${BOOTSTRAPS}"
    echo "focal_species=${KSRATES_ONLY}"
    echo "primary_haplotypes=P. zippeliana A1+B1; P. spinulosa A1+B1; P. padus A+B"
    echo "reference_species=P. campanulata,P. persica,G. trifoliata"
} > "${KS_ROOT}/ksrates/software_and_parameters.txt"

all_complete=1
for focal in pzi psp ppa; do
    [[ -e "${RUN_ROOT}/${focal}/COMPLETE" ]] || all_complete=0
done
if [[ "${all_complete}" == "1" ]]; then
    touch "${KS_ROOT}/ksrates/COMPLETE"
    echo "[OK] All three ksrates analyses completed: ${RUN_ROOT}"
else
    echo "[OK] Requested ksrates analyses completed; other focal runs remain pending: ${RUN_ROOT}"
fi
