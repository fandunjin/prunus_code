#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FIG_ROOT="${FIG_ROOT:-${PAN_ROOT}/figure1_tandem}"
PYTHON="${PYTHON:-${CONDA_ROOT}/envs/haphic/bin/python}"
SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_tandem_loci.py"
LOCI="${FIG_ROOT}/loci"
TANDEM_ROOT="${PAN_ROOT}/followup/tandem"
MANIFEST="${LOCI}/pzi_tandem_flank_comparisons_manifest.tsv"

for path in "${PYTHON}" "${SCRIPT}" "${LOCI}/pzi_tandem_array_candidates.tsv" "${TANDEM_ROOT}"; do
    [ -e "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
for target in pziA1 pziA2 pziB1 pziB2 pspA1 pspA2 pspB1 pspB2 ppaA ppaB pcam ppe; do
    [ -s "${LOCI}/paf/all_pzi_array_flanks_to_${target}.paf" ] || { echo "ERROR: missing PAF for ${target}" >&2; exit 1; }
done

printf 'query_target\ttarget\tpaf\n' > "${MANIFEST}"
for query in pziA1 pziA2; do
    if [ "${query}" = pziA1 ]; then counterpart=pziA2; else counterpart=pziA1; fi
    printf '%s\t%s\t%s\n' "${query}" "${counterpart}" "${LOCI}/paf/all_pzi_array_flanks_to_${counterpart}.paf" >> "${MANIFEST}"
    for target in pspA1 pspA2; do
        printf '%s\t%s\t%s\n' "${query}" "${target}" "${LOCI}/paf/all_pzi_array_flanks_to_${target}.paf" >> "${MANIFEST}"
    done
    for target in ppaA ppaB pcam ppe; do
        printf '%s\t%s\t%s\n' "${query}" "${target}" "${LOCI}/paf/all_pzi_array_flanks_to_${target}.paf" >> "${MANIFEST}"
    done
done
for query in pziB1 pziB2; do
    if [ "${query}" = pziB1 ]; then counterpart=pziB2; else counterpart=pziB1; fi
    printf '%s\t%s\t%s\n' "${query}" "${counterpart}" "${LOCI}/paf/all_pzi_array_flanks_to_${counterpart}.paf" >> "${MANIFEST}"
    for target in pspB1 pspB2; do
        printf '%s\t%s\t%s\n' "${query}" "${target}" "${LOCI}/paf/all_pzi_array_flanks_to_${target}.paf" >> "${MANIFEST}"
    done
    for target in ppaA ppaB pcam ppe; do
        printf '%s\t%s\t%s\n' "${query}" "${target}" "${LOCI}/paf/all_pzi_array_flanks_to_${target}.paf" >> "${MANIFEST}"
    done
done

"${PYTHON}" "${SCRIPT}" summarize-flanks \
    --candidates "${LOCI}/pzi_tandem_array_candidates.tsv" \
    --comparison-manifest "${MANIFEST}" \
    --tandem-root "${TANDEM_ROOT}" \
    --output-dir "${LOCI}" \
    --min-trf-score 50 \
    --min-mapq 20 \
    --min-query-coverage 0.80 \
    --min-identity 0.70 \
    --contracted-ratio 0.50 \
    --expanded-ratio 1.50

echo "=== Tandem-array locus summary completed ==="
