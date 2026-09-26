#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# QUERY_TARGET=pziB1 jsub < jobs/haploid_20a_prunus_pan_te_orthosite.sh

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/prunus_pan_te}"
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/prunus_pan_te_orthosite.py}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
DECISIONS="${DECISIONS:-${FOLLOWUP_ROOT}/component_audit/component_decisions.tsv}"
INTACT_LTR_TABLE="${INTACT_LTR_TABLE:-${OUT_ROOT}/postprocess_v2/te_intact_ltr.tsv}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MAP_ENV="${MAP_ENV:-minimap}"
THREADS="${THREADS:-8}"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-}"
MINIMAP_PRESET="${MINIMAP_PRESET:-asm20}"
MAX_COMPONENTS="${MAX_COMPONENTS:-10}"
MAX_LOCI_PER_COMPONENT="${MAX_LOCI_PER_COMPONENT:-200}"
FLANK_BP="${FLANK_BP:-1500}"
FORCE_RERUN="${FORCE_RERUN:-0}"

case "${QUERY_TARGET}" in
    pziA1)
        # Keep all PSP/PPA haplotypes in the comparison set.  A/B labels are
        # species-specific working labels and must not be assumed homologous
        # before the cross-species subgenome test.
        COMPARISONS=(
            "pziA2:within_pzi"
            "pspA1:psp" "pspA2:psp" "pspB1:psp" "pspB2:psp"
            "ppaA:outgroup" "ppaB:outgroup"
            "pcam:outgroup" "ppe:outgroup"
        )
        ;;
    pziB1)
        COMPARISONS=(
            "pziB2:within_pzi"
            "pspA1:psp" "pspA2:psp" "pspB1:psp" "pspB2:psp"
            "ppaA:outgroup" "ppaB:outgroup"
            "pcam:outgroup" "ppe:outgroup"
        )
        ;;
    *)
        echo "ERROR: QUERY_TARGET must be pziA1 or pziB1" >&2
        exit 1
        ;;
esac

for path in "${SCRIPT}" "${TARGET_CONFIG}" "${DECISIONS}" "${INTACT_LTR_TABLE}" "${CONDA_SH}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${MAP_ENV}"
for exe in minimap2 samtools python; do
    command -v "${exe}" >/dev/null 2>&1 || { echo "ERROR: ${exe} was not found" >&2; exit 1; }
done
if [ -z "${BEDTOOLS_BIN}" ]; then
    BEDTOOLS_BIN="$(command -v bedtools || true)"
elif [[ "${BEDTOOLS_BIN}" != */* ]]; then
    BEDTOOLS_BIN="$(command -v "${BEDTOOLS_BIN}" || true)"
fi
[ -n "${BEDTOOLS_BIN}" ] && [ -x "${BEDTOOLS_BIN}" ] || {
    echo "ERROR: bedtools was not found; set BEDTOOLS_BIN to an executable path" >&2
    exit 1
}

query_genome="${OUT_ROOT}/mask/${QUERY_TARGET}/${QUERY_TARGET}.fa"
query_rmout="${query_genome}.out"
outdir="${FOLLOWUP_ROOT}/orthosite/${QUERY_TARGET}"
mkdir -p "${outdir}/paf"
for path in "${query_genome}" "${query_rmout}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
[ -s "${query_genome}.fai" ] || samtools faidx "${query_genome}"

python "${SCRIPT}" prepare \
    --query-target "${QUERY_TARGET}" \
    --rmout "${query_rmout}" \
    --fai "${query_genome}.fai" \
    --decisions "${DECISIONS}" \
    --intact-ltr-table "${INTACT_LTR_TABLE}" \
    --output-dir "${outdir}" \
    --max-components "${MAX_COMPONENTS}" \
    --max-loci-per-component "${MAX_LOCI_PER_COMPONENT}" \
    --flank-bp "${FLANK_BP}"

"${BEDTOOLS_BIN}" getfasta \
    -fi "${query_genome}" \
    -bed "${outdir}/query_flanks.bed" \
    -nameOnly \
    -fo "${outdir}/query_flanks.fa"

manifest="${outdir}/comparison_targets.tsv"
printf "target\trole\tgenome\trepeatmasker_out\tpaf\tminimap_preset\n" > "${manifest}"
for item in "${COMPARISONS[@]}"; do
    target="${item%%:*}"
    role="${item##*:}"
    target_genome="${OUT_ROOT}/mask/${target}/${target}.fa"
    target_rmout="${target_genome}.out"
    paf="${outdir}/paf/${QUERY_TARGET}_to_${target}.flanks.paf"
    for path in "${target_genome}" "${target_rmout}"; do
        [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
    done
    if [ "${FORCE_RERUN}" = "1" ] || [ ! -s "${paf}" ]; then
        minimap2 -x "${MINIMAP_PRESET}" -N 5 --secondary=yes -t "${THREADS}" \
            "${target_genome}" "${outdir}/query_flanks.fa" > "${paf}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${target}" "${role}" "${target_genome}" "${target_rmout}" "${paf}" \
        "${MINIMAP_PRESET}" >> "${manifest}"
done

python "${SCRIPT}" summarize \
    --loci "${outdir}/query_loci.tsv" \
    --target-manifest "${manifest}" \
    --decisions "${DECISIONS}" \
    --output-dir "${outdir}"

echo "=== Orthologous insertion-site analysis completed: ${QUERY_TARGET} ==="
