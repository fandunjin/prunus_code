#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# te_sololtr_panel10_submit_all.sh
# Preflight + submit 10 soloLTRseeker jobs (panel-10: pziA1/B1, pspA1/B1, ppaA/B, pcam, ppe, gtr, rch).
# Run on the HPC login node:  bash jobs/te_sololtr_panel10_submit_all.sh
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
SOLOLTR_ENV="${SOLOLTR_ENV:-soloLTR}"

TARGETS=(pziA1 pziB1 pspA1 pspB1 ppaA ppaB pcam ppe gtr rch)

die() { echo "[ERROR] $*" >&2; exit 2; }

echo "=== conda env preflight (${SOLOLTR_ENV}) ==="
[[ -s "${CONDA_SH}" ]] || die "missing ${CONDA_SH}"
source "${CONDA_SH}"
conda activate "${SOLOLTR_ENV}"
for exe in bedtools blastn cd-hit-est water TEsorter; do
    command -v "${exe}" >/dev/null 2>&1 || die "${exe} not found in env ${SOLOLTR_ENV}"
done
SOLOLTR_BIN="${SOLOLTR_BIN:-}"
if [ -z "${SOLOLTR_BIN}" ]; then
    for candidate in \
        "${ROOT_DIR}/software/soloLTRseeker/modules/soloLTRseeker" \
        "${ROOT_DIR}/software/soloLTRseeker/soloLTRseeker" \
        "${ROOT_DIR}/software/soloLTRseeker/bin/soloLTRseeker" \
        "${ROOT_DIR}/../software/soloLTRseeker-main/modules/soloLTRseeker" \
        "${ROOT_DIR}/../software/soloLTRseeker-main/soloLTRseeker" \
        "${ROOT_DIR}/../software/soloLTRseeker-main/bin/soloLTRseeker" \
        "${ROOT_DIR}/../software/soloLTRseeker/modules/soloLTRseeker"; do
        if [ -s "${candidate}" ]; then SOLOLTR_BIN="${candidate}"; break; fi
    done
fi
if [ -z "${SOLOLTR_BIN}" ] && [ -d "${ROOT_DIR}/../software" ]; then
    SOLOLTR_BIN="$(find "${ROOT_DIR}/../software" -maxdepth 4 -type f -name 'soloLTRseeker' 2>/dev/null | head -n 1)"
fi
if [ -z "${SOLOLTR_BIN}" ]; then
    SOLOLTR_BIN="$(command -v soloLTRseeker || true)"
fi
if [ -n "${SOLOLTR_BIN}" ] && [ -s "${SOLOLTR_BIN}" ]; then
    echo "  soloLTRseeker: ${SOLOLTR_BIN}"
else
    die "soloLTRseeker not found (set SOLOLTR_BIN or install under software/soloLTRseeker)"
fi

echo "=== input preflight ==="
missing=0
for t in "${TARGETS[@]}"; do
    g=""
    case "${t}" in
        pziA1|pziB1)
            h="${t#pzi}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/pzi/edta/hap${h}/pzi_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        pspA1|pspB1)
            h="${t#psp}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/psp/edta/hap${h}/psp_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        ppaA|ppaB)
            h="${t#ppa}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/ppa/edta/hap${h}/ppa_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        pcam|ppe|gtr|rch)
            for base in "${ROOT_DIR}/TE" "${ROOT_DIR}/0903/TE"; do
                cand="${base}/${t}/${t}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
    esac
    if [ -n "${g}" ] && [ -s "${g}" ] && [ -s "${g}.mod.EDTA.intact.gff3" ]; then
        printf '  %-6s OK   %s\n' "${t}" "${g}"
    else
        printf '  %-6s MISSING (genome/intact GFF)\n' "${t}"
        missing=$((missing + 1))
    fi
done
[[ "${missing}" -eq 0 ]] || die "${missing} input(s) missing; fix paths first"

if [[ "${PREFLIGHT_ONLY:-0}" == "1" ]]; then
    echo "=== PREFLIGHT_ONLY=1 : inputs and environment OK, submission skipped ==="
    exit 0
fi

cd "${ROOT_DIR}"
for t in "${TARGETS[@]}"; do
    export TARGET="${t}"
    echo ">>> submit ${t}"
    jsub < jobs/te_sololtr_panel10_run.sh
done
echo "=== 10 soloLTRseeker jobs submitted ==="
