#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Neutral copy1-copy2 and independent haplotype synteny.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_internal_synteny}"
THREADS="${THREADS:-12}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

COPY1_DIR="${COPY1_DIR:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison/inputs/pzi_copy1}"
COPY2_GENOME="${COPY2_GENOME:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy2.hap_provenance.chr16.fa}"
HAP_INPUT_ROOT="${HAP_INPUT_ROOT:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs}"
COPY_MAP="${COPY_MAP:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv}"
GTR_PROFILES="${GTR_PROFILES:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison/model_results/chromosome_profiles.tsv}"

# Prefer the current EDTA-masked Helixer annotation.  The fallback keeps this
# driver usable with the earlier unmasked copy2 annotation if that result is
# the only one already present on the cluster.
if [[ -z "${COPY2_GFF:-}" ]]; then
    for candidate in \
        "${ROOT_DIR}/annotation/genes/copy2_helixer_edta_masked_optional/pzi_copy2.edta_masked.helixer.gff3" \
        "${ROOT_DIR}/annotation/genes/copy2_helixer_optional/pzi_copy2.helixer.gff3"; do
        if [[ -s "${candidate}" ]]; then COPY2_GFF="${candidate}"; break; fi
    done
fi
COPY2_GFF="${COPY2_GFF:-${ROOT_DIR}/annotation/genes/copy2_helixer_edta_masked_optional/pzi_copy2.edta_masked.helixer.gff3}"

# This summary is produced by final_homolog_syri_tidk.sh.  Accept the current
# final_assembly location and the packaged copy for portability.
if [[ -z "${SAME_CHR_ALIGNMENT:-}" ]]; then
    for candidate in \
        "${ROOT_DIR}/final_assembly/homolog_syri_tidk/pzi.copy1_vs_copy2.same_chr_alignment.tsv" \
        "${ROOT_DIR}/final_assembly/pzi.copy1_vs_copy2.same_chr_alignment.tsv" \
        "${ROOT_DIR}/final_assembly/haplotype_provenance/pzi.copy1_vs_copy2.same_chr_alignment.tsv"; do
        if [[ -s "${candidate}" ]]; then SAME_CHR_ALIGNMENT="${candidate}"; break; fi
    done
fi
SAME_CHR_ALIGNMENT="${SAME_CHR_ALIGNMENT:-${ROOT_DIR}/final_assembly/homolog_syri_tidk/pzi.copy1_vs_copy2.same_chr_alignment.tsv}"

INPUT_ROOT="${OUT_ROOT}/inputs"
PAIR_ROOT="${OUT_ROOT}/pairs"
SUMMARY_ROOT="${OUT_ROOT}/summary"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing input: $1"; }

for f in "${COPY1_DIR}/pzi_copy1.bed" "${COPY1_DIR}/pzi_copy1.pep" "${COPY1_DIR}/pzi_copy1.cds" \
         "${COPY2_GFF}" "${COPY2_GENOME}" "${COPY_MAP}" "${SAME_CHR_ALIGNMENT}" "${GTR_PROFILES}"; do req "${f}"; done
for src in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2; do
    req "${HAP_INPUT_ROOT}/${src}/${src}.bed"
    req "${HAP_INPUT_ROOT}/${src}/${src}.pep"
    req "${HAP_INPUT_ROOT}/${src}/${src}.cds"
done
if [[ -e "${OUT_ROOT}/COMPLETE" && "${FORCE}" != "1" ]]; then
    echo "[OK] Existing neutral internal synteny: ${OUT_ROOT}"
    exit 0
fi

mkdir -p "${INPUT_ROOT}" "${PAIR_ROOT}" "${SUMMARY_ROOT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_neutral_internal_pzi_inputs.py" \
    --copy1-dir "${COPY1_DIR}" \
    --copy2-gff "${COPY2_GFF}" \
    --copy2-genome "${COPY2_GENOME}" \
    --hap hap01 "${HAP_INPUT_ROOT}/pzi_hapA1" \
    --hap hap02 "${HAP_INPUT_ROOT}/pzi_hapA2" \
    --hap hap03 "${HAP_INPUT_ROOT}/pzi_hapB1" \
    --hap hap04 "${HAP_INPUT_ROOT}/pzi_hapB2" \
    --outdir "${INPUT_ROOT}"

for label in copy1 copy2 hap01 hap02 hap03 hap04; do
    req "${INPUT_ROOT}/${label}/${label}.bed"
    req "${INPUT_ROOT}/${label}/${label}.pep"
    req "${INPUT_ROOT}/${label}/${label}.cds"
done

labels=(copy1 copy2 hap01 hap02 hap03 hap04)
plan="${SUMMARY_ROOT}/pair_plan.tsv"
{
    printf 'comparison\tquery_label\tref_label\tcomparison_type\n'
    for ((i=0; i<${#labels[@]}; i++)); do
        for ((j=i+1; j<${#labels[@]}; j++)); do
            q="${labels[$i]}"; r="${labels[$j]}"
            type="independent_haplotype"
            [[ "${q}_${r}" == "copy1_copy2" ]] && type="copy1_copy2_internal"
            [[ "${q}" == copy1 || "${q}" == copy2 ]] && [[ "${r}" == hap* ]] && type="copy_to_haplotype"
            printf '%s\t%s\t%s\t%s\n' "${q}_vs_${r}" "${q}" "${r}" "${type}"
        done
    done
} > "${plan}"

command -v diamond >/dev/null 2>&1 || die "diamond is unavailable on PATH"
"${PYTHON_BIN}" -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is unavailable in ${PYTHON_BIN}"

while IFS=$'\t' read -r comparison query ref comparison_type; do
    [[ "${comparison}" == "comparison" ]] && continue
    work="${PAIR_ROOT}/${comparison}"
    mkdir -p "${work}"
    for label in "${query}" "${ref}"; do
        for ext in bed pep cds; do
            ln -sfn "${INPUT_ROOT}/${label}/${label}.${ext}" "${work}/${label}.${ext}"
        done
    done
    pushd "${work}" >/dev/null
    raw="${query}.${ref}.anchors"
    screened="${query}.${ref}.anchors.new"
    if [[ "${FORCE}" == "1" ]]; then rm -f "${raw}" "${screened}" "${screened}.simple" "${query}.${ref}.anchors.simple"; fi
    if [[ ! -s "${raw}" ]]; then
        "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
            --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" \
            --min_size="${MIN_ANCHORS}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    fi
    req "${work}/${raw}"
    if [[ ! -s "${screened}" ]]; then
        "${PYTHON_BIN}" -m jcvi.compara.synteny screen \
            --qbed="${query}.bed" --sbed="${ref}.bed" --minspan="${MINSPAN}" \
            --simple "${raw}" "${screened}"
    fi
    req "${work}/${screened}"
    popd >/dev/null
done < "${plan}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/summarize_neutral_internal_synteny.py" \
    --pair-plan "${plan}" \
    --pair-root "${PAIR_ROOT}" \
    --input-root "${INPUT_ROOT}" \
    --copy-map "${COPY_MAP}" \
    --same-chr-alignment "${SAME_CHR_ALIGNMENT}" \
    --gtr-profiles "${GTR_PROFILES}" \
    --outdir "${SUMMARY_ROOT}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/integrate_neutral_internal_model_evidence.py" \
    --model-comparison "${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison/model_results/model_comparison.tsv" \
    --k97-partitions "${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_copy1_gtr_model_comparison/model_results/k9_plus_7_top_partitions.tsv" \
    --resolution "${SUMMARY_ROOT}/copy1_copy2_gtr_resolution.tsv" \
    --hap-pairs "${SUMMARY_ROOT}/independent_haplotype_pairs.tsv" \
    --outdir "${SUMMARY_ROOT}/model_integration"

{
    echo "python=$(${PYTHON_BIN} --version 2>&1)"
    echo "jcvi=$(${PYTHON_BIN} -c 'import jcvi; print(getattr(jcvi, "__version__", "unknown"))' 2>/dev/null || echo unavailable)"
    echo "diamond=$(diamond version 2>&1 | head -n 1)"
    echo "cscore=${CSCORE}"
    echo "min_size=${MIN_ANCHORS}"
    echo "minspan=${MINSPAN}"
    echo "datasets=copy1,copy2,hap01,hap02,hap03,hap04"
    echo "interpretation=neutral internal synteny; source hap names are metadata only"
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Neutral internal synteny completed: ${OUT_ROOT}"
