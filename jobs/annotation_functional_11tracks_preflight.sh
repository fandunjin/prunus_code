#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Preflight check for the 11-track Swiss-Prot annotation jobs.
# Run on the HPC login node BEFORE submitting:
#   cd ${PROJECT_ROOT}/annotation_functional_11tracks_preflight.sh
# Then:
#   jsub < jobs/annotation_functional_11tracks.sh
#   jsub < jobs/annotation_functional_11tracks_metadata.sh
#
# This script only reads state; it never modifies inputs or outputs.

set -uo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
FUNCTIONAL_ENV="${FUNCTIONAL_ENV:-orthofinder}"
DIAMOND_BIN="${DIAMOND_BIN:-diamond}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SWISSPROT_DB="${SWISSPROT_DB:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
SWISSPROT_DAT="${SWISSPROT_DAT:-${HOME}/CX/database/SwissProt/uniprot_sprot.dat.gz}"
FROZEN_PEP="${FROZEN_PEP:-${ROOT_DIR}/annotation/functional/inputs_11tracks}"
FUNCTIONAL_ROOT="${FUNCTIONAL_ROOT:-${ROOT_DIR}/annotation/functional}"
WORK_DMND="${FUNCTIONAL_ROOT}/work/swissprot.dmnd"
META_SCRIPT="${ROOT_DIR}/scripts/extract_swissprot_metadata.py"
JOB1="${ROOT_DIR}/jobs/annotation_functional_11tracks.sh"
JOB2="${ROOT_DIR}/jobs/annotation_functional_11tracks_metadata.sh"
TRACKS="${SPECIES_SET:-pcam pav pspe pmal par pmu pdu psal ppe gtr rch}"
LOG_DIR="${ROOT_DIR}/logs"

PASS=0; FAIL=0; WARN=0
ok()   { printf '  [ OK ]  %s\n' "$*"; PASS=$((PASS+1)); }
fail() { printf '  [FAIL]  %s\n' "$*"; FAIL=$((FAIL+1)); }
warn() { printf '  [WARN]  %s\n' "$*"; WARN=$((WARN+1)); }
info() { printf '  [info]  %s\n' "$*"; }

echo "=================================================================="
echo "Preflight: 11-track Swiss-Prot annotation"
echo "ROOT_DIR=${ROOT_DIR}"
echo "=================================================================="

echo
echo "1) Software and conda environment"
if [[ -s "${CONDA_SH}" ]]; then ok "conda init script: ${CONDA_SH}"; else fail "conda init script missing: ${CONDA_SH}"; fi
if [[ -s "${CONDA_SH}" ]] && source "${CONDA_SH}" >/dev/null 2>&1 && conda activate "${FUNCTIONAL_ENV}" >/dev/null 2>&1; then
    ok "conda env '${FUNCTIONAL_ENV}' activates"
else
    fail "cannot activate conda env '${FUNCTIONAL_ENV}' (check CONDA_SH / env name)"
fi
if command -v "${DIAMOND_BIN}" >/dev/null 2>&1; then
    ok "${DIAMOND_BIN}: $(command -v "${DIAMOND_BIN}")  [$( "${DIAMOND_BIN}" version 2>&1 | head -n 1 )]"
else
    fail "${DIAMOND_BIN} not on PATH inside env '${FUNCTIONAL_ENV}'"
fi
if command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    ok "${PYTHON_BIN}: $(command -v "${PYTHON_BIN}")  [$( "${PYTHON_BIN}" --version 2>&1 )]"
else
    fail "${PYTHON_BIN} not on PATH inside env '${FUNCTIONAL_ENV}'"
fi

echo
echo "2) Job scripts and metadata parser"
for f in "${JOB1}" "${JOB2}" "${META_SCRIPT}"; do
    if [[ -s "${f}" ]]; then ok "$(basename "${f}")  ($(stat -c '%s bytes' "${f}" 2>/dev/null || echo size?))"; else fail "missing: ${f}"; fi
done

echo
echo "3) Swiss-Prot databases"
if [[ -s "${WORK_DMND}" ]]; then
    ok "prebuilt DIAMOND db: ${WORK_DMND}  ($(stat -c '%s bytes' "${WORK_DMND}" 2>/dev/null || echo size?))"
else
    warn "prebuilt db not found: ${WORK_DMND} (job will build it from SWISSPROT_DB, needs time/disk)"
fi
if [[ -s "${SWISSPROT_DB}" ]]; then
    ok "Swiss-Prot FASTA: ${SWISSPROT_DB}  ($(stat -c '%s bytes' "${SWISSPROT_DB}" 2>/dev/null || echo size?))"
elif [[ "${SWISSPROT_DB}" == *.dmnd && -s "${SWISSPROT_DB}" ]]; then
    ok "Swiss-Prot dmnd: ${SWISSPROT_DB}"
else
    fail "Swiss-Prot database missing: ${SWISSPROT_DB}"
fi
if [[ -s "${SWISSPROT_DAT}" ]]; then
    ok "Swiss-Prot flat file: ${SWISSPROT_DAT}  ($(stat -c '%s bytes' "${SWISSPROT_DAT}" 2>/dev/null || echo size?))"
else
    fail "Swiss-Prot flat file missing: ${SWISSPROT_DAT} (needed by stage 2)"
fi
if [[ -s "${SWISSPROT_DB}" && -s "${SWISSPROT_DAT}" ]]; then
    warn "confirm uniprot_sprot.fasta and uniprot_sprot.dat.gz are the SAME release (compare download date/md5)"
fi

echo
echo "4) Frozen protein inputs for the 11 tracks (species-level tracks only)"
# SHA256 of the current-revision frozen inputs (identical in both desktop mirrors).
# Any mismatch means the file may come from a superseded revision.
declare -A EXPECTED_SHA=(
    [pcam]="1219191708e9c700feff2ebd433c27db084774431b26da4b45050310944168d1"
    [pav]="8f94b88eb243ae115e9d14acf4577a529eaaa326ceaf3316a2f558d5f290079d"
    [pspe]="9f1f44ba783ebd9234f90d44d1d62f6f4fc15d2742cc909236aad8b686f72f5e"
    [pmal]="d6e33f88ff1167883f83a93a697807983df5e6dbee2b77902e831652e3341669"
    [par]="b0f08a1a7d76e4b8bf71915f041e079cf1bb3c78b84e433a9a720740b251f577"
    [pmu]="f10e0b31cb1fd623c699bd231f496906d41a940ada021d968df77881e79db7b4"
    [pdu]="bc8411f286aaa8f5f34be2f1cad8f8ea63ae7321cbb53c9566053386330abe97"
    [psal]="ef0d330caff466b4955a284ffc096bda93b84ca2f36e4218ff4852f5bdeaefbe"
    [ppe]="73716618ca73d02960a8e0122160b572e2a2edeb963823cf0377e16e93bd33f6"
    [gtr]="8b7dfc5fe2c93920b18ce914f3078b3a641987a862a3b6db9b6e22ea7edced0a"
    [rch]="a39468ae92b7d4be8293dac2dd91d6f78b98ebb22d31a7e5d52b1ed91dd2566e"
)
SHA_BIN="$(command -v sha256sum || true)"
[[ -n "${SHA_BIN}" ]] || warn "sha256sum not available; hash verification skipped"
PEP_CANDIDATES=(
    "${FROZEN_PEP}"
    "${ROOT_DIR}/downstream_analysis/phylogenetics/01_frozen_inputs/main/pep"
    "${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs/main/pep"
    "${ROOT_DIR}/analysisphylogenetics/01_frozen_inputs/main/pep"
)
MISSING_PEP=0
for t in ${TRACKS}; do
    pep=""
    for d in "${PEP_CANDIDATES[@]}"; do
        [[ -s "${d}/${t}.fa" ]] && pep="${d}/${t}.fa" && break
    done
    if [[ -n "${pep}" ]]; then
        n=$(grep -c '^>' "${pep}" 2>/dev/null || echo "?")
        size=$(stat -c '%s bytes' "${pep}" 2>/dev/null || echo "size?")
        if [[ -n "${SHA_BIN}" ]]; then
            got=$("${SHA_BIN}" "${pep}" | awk '{print $1}')
            if [[ "${got}" == "${EXPECTED_SHA[$t]}" ]]; then
                ok "${t}: ${pep}  (${n} proteins, ${size}) sha256 OK"
            else
                warn "${t}: sha256 MISMATCH - got ${got:0:12}..., expected ${EXPECTED_SHA[$t]:0:12}...; file may be from a superseded revision"
            fi
        else
            ok "${t}: ${pep}  (${n} proteins, ${size})"
        fi
    else
        fail "${t}: protein FASTA not found in any of: ${PEP_CANDIDATES[*]}"
        MISSING_PEP=1
    fi
done
if [[ "${MISSING_PEP}" == 1 ]]; then
    info "discovery hints:"
    info "  ls -d ${ROOT_DIR}/*/ | head -40"
    info "  find ${ROOT_DIR} -maxdepth 7 -type d -name '01_frozen_inputs' 2>/dev/null"
    info "  find ${ROOT_DIR} -maxdepth 8 -type f -name 'pcam.fa' 2>/dev/null | head"
    if command -v find >/dev/null 2>&1; then
        FOUND="$(find "${ROOT_DIR}" -maxdepth 7 -type d -name '01_frozen_inputs' 2>/dev/null | head -n 5)"
        if [[ -n "${FOUND}" ]]; then
            info "find results for 01_frozen_inputs:"
            while IFS= read -r line; do info "  ${line}"; done <<< "${FOUND}"
        fi
    fi
fi

echo
echo "5) Output locations (existing results make jobs SKIP unless FORCE=1)"
if mkdir -p "${FUNCTIONAL_ROOT}" "${LOG_DIR}" 2>/dev/null && [[ -w "${FUNCTIONAL_ROOT}" && -w "${LOG_DIR}" ]]; then
    ok "writable: ${FUNCTIONAL_ROOT} and ${LOG_DIR}"
else
    fail "not writable: ${FUNCTIONAL_ROOT} or ${LOG_DIR}"
fi
for t in ${TRACKS}; do
    if [[ -s "${FUNCTIONAL_ROOT}/${t}/swissprot_functional_annotation.tsv" ]]; then
        warn "${t}: swissprot_functional_annotation.tsv already exists (FORCE=1 to rerun)"
    elif [[ -s "${FUNCTIONAL_ROOT}/${t}/swissprot.diamond.tsv" ]]; then
        info "${t}: stage-1 output exists, stage 2 pending"
    else
        info "${t}: no existing output"
    fi
done
for t in pzi psp ppa; do
    if [[ -s "${FUNCTIONAL_ROOT}/${t}/swissprot_functional_annotation.tsv" ]]; then
        info "guard: ${t} annotation present and will NOT be touched by these jobs"
    fi
done

echo
echo "6) Disk space"
if command -v df >/dev/null 2>&1; then
    df -h "${ROOT_DIR}" | tail -n 1 | awk '{printf "  [info]  %s available on %s\n", $4, $6}'
fi

echo
echo "=================================================================="
printf "Summary: %d OK, %d WARN, %d FAIL\n" "${PASS}" "${WARN}" "${FAIL}"
if [[ "${FAIL}" -gt 0 ]]; then
    echo "STATUS: NOT READY - fix the [FAIL] items before submitting."
    exit 1
fi
echo "STATUS: READY - submit with:"
echo "  cd ${ROOT_DIR}"
echo "  jsub < jobs/annotation_functional_11tracks.sh"
echo "  # after stage 1 finishes (check logs/annotation_functional_11tracks_*.log):"
echo "  jsub < jobs/annotation_functional_11tracks_metadata.sh"
echo "=================================================================="
