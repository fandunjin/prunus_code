#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# te_sololtr_panel10_run.sh
# soloLTRseeker for one panel-10 target (pziA1/B1, pspA1/B1, ppaA/B, pcam, ppe, gtr, rch).
# Submit from the repo root with TARGET exported, e.g.
#   TARGET=gtr jsub < jobs/te_sololtr_panel10_run.sh
# Optional: CLEAN=1 removes an existing output dir first.
SOLOLTR_DIR="${SOLOLTR_DIR:-${ROOT_DIR}/software/soloLTRseeker/modules}"
SOLOLTR_BIN="${SOLOLTR_BIN:-}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SOLOLTR_ENV="${SOLOLTR_ENV:-soloLTR}"
MIN_LTR_LENGTH="${MIN_LTR_LENGTH:-100}"
MAX_LTR_PAIR_LENGTH_DIFF="${MAX_LTR_PAIR_LENGTH_DIFF:-50}"

die() { echo "[ERROR] $*" >&2; exit 2; }

case "${TARGET}" in
    pziA1|pziB1|pspA1|pspB1|ppaA|ppaB|pcam|ppe|gtr|rch) ;;
    *) die "unsupported TARGET='${TARGET}' (expected one of the panel-10 targets)" ;;
esac

resolve_inputs() {
    local g="" cand="" base=""
    case "${TARGET}" in
        pziA1|pziB1)
            local h="${TARGET#pzi}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/pzi/edta/hap${h}/pzi_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        pspA1|pspB1)
            local h="${TARGET#psp}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/psp/edta/hap${h}/psp_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        ppaA|ppaB)
            local h="${TARGET#ppa}"
            for base in "${ROOT_DIR}/haploid_repartition" "${ROOT_DIR}/0903/haploid_repartition"; do
                cand="${base}/inputs/ppa/edta/hap${h}/ppa_hap${h}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
        pcam|ppe|gtr|rch)
            for base in "${ROOT_DIR}/TE" "${ROOT_DIR}/0903/TE"; do
                cand="${base}/${TARGET}/${TARGET}.fa"
                if [ -s "${cand}" ]; then g="${cand}"; break; fi
            done
            ;;
    esac
    [ -n "${g}" ] && [ -s "${g}" ] || die "genome FASTA not found for ${TARGET}"
    local ig="${g}.mod.EDTA.intact.gff3"
    [ -s "${ig}" ] || die "intact GFF3 not found: ${ig}"
    local gfm="${g%.fa}.gffmatch.fa"
    if [ -n "${GENOME_OVERRIDE}" ]; then
        [ -s "${GENOME_OVERRIDE}" ] || die "GENOME_OVERRIDE not found: ${GENOME_OVERRIDE}"
        GENOME="${GENOME_OVERRIDE}"
    elif [ -s "${gfm}" ]; then
        GENOME="${gfm}"
        echo "[INFO] using GFF-compatible genome FASTA: ${gfm}"
    else
        GENOME="${g}"
    fi
    INTACT_GFF="${INTACT_GFF_OVERRIDE:-${ig}}"
    [ -s "${INTACT_GFF}" ] || die "intact GFF3 not found: ${INTACT_GFF}"
}

[[ -s "${CONDA_SH}" ]] || die "missing ${CONDA_SH}"
resolve_inputs
source "${CONDA_SH}"
conda activate "${SOLOLTR_ENV}"

if [ -z "${SOLOLTR_BIN}" ]; then
    for candidate in \
        "${SOLOLTR_DIR}/soloLTRseeker" \
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
[ -n "${SOLOLTR_BIN}" ] && [ -s "${SOLOLTR_BIN}" ] || die "soloLTRseeker not found (set SOLOLTR_BIN)"
[ -x "${SOLOLTR_BIN}" ] || die "soloLTRseeker is not executable: ${SOLOLTR_BIN} (run: chmod u+x ${SOLOLTR_BIN})"
for exe in bedtools blastn cd-hit-est water TEsorter; do
    command -v "${exe}" >/dev/null 2>&1 || die "${exe} not found in env ${SOLOLTR_ENV}"
done

outdir="${OUT_ROOT}/sololtr/${TARGET}"
if [[ "${CLEAN}" == "1" && -d "${outdir}" ]]; then
    echo "[CLEAN] removing ${outdir}"
    rm -rf "${outdir}"
fi
mkdir -p "${outdir}" "${ROOT_DIR}/logs"
cd "${outdir}"

echo "=== soloLTRseeker: ${TARGET} ==="
echo "Genome      : ${GENOME}"
echo "Intact GFF3 : ${INTACT_GFF}"
echo "Output      : ${outdir}"
echo "Cores       : ${THREADS}"

sololtr_genome="${outdir}/${TARGET}.fa"
if [ ! -s "${sololtr_genome}" ]; then
    tmp_genome="${sololtr_genome}.tmp.$$"
    awk '/^>/{sub(/[[:space:]].*$/, "")}1' "${GENOME}" > "${tmp_genome}"
    [ -s "${tmp_genome}" ] || die "failed to create ${tmp_genome}"
    mv "${tmp_genome}" "${sololtr_genome}"
fi

if ! find . -maxdepth 1 -type f -name '*_soloLTR.gff3' -size +0c | grep -q .; then
    "${SOLOLTR_BIN}" \
        -n "${MIN_LTR_LENGTH}" \
        -b "${MAX_LTR_PAIR_LENGTH_DIFF}" \
        -t \
        -c 0.95 \
        -l 0.90 \
        -s 0.30 \
        -o 0.99 \
        -q 200 \
        -p 6 \
        -m 1 \
        -r "${THREADS}" \
        -u T \
        "${INTACT_GFF}" "${sololtr_genome}"
fi

solo_file="$(find . -type f -name '*_soloLTR.gff3' -size +0c -printf '%T@\t%p\n' | sort -nr | head -n 1 | cut -f2-)"
length_file="$(find . -type f -name '*_LTRRT_lengths.txt' -size +0c -printf '%T@\t%p\n' | sort -nr | head -n 1 | cut -f2-)"
[ -n "${solo_file}" ] || die "no non-empty soloLTR GFF found under ${outdir}"
[ -n "${length_file}" ] || die "no non-empty LTRRT length table found under ${outdir}"
ln -sfn "$(readlink -f "${solo_file}")" "${TARGET}.soloLTR.gff3"
ln -sfn "$(readlink -f "${length_file}")" "${TARGET}.LTRRT_lengths.txt"

printf 'target\tgenome\tintact_gff\tsolo_gff\tltrrt_lengths\tmin_ltr_length\tmax_pair_length_difference\n' > "${TARGET}.sololtr_manifest.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${TARGET}" "${GENOME}" "${INTACT_GFF}" "$(readlink -f "${TARGET}.soloLTR.gff3")" \
    "$(readlink -f "${TARGET}.LTRRT_lengths.txt")" "${MIN_LTR_LENGTH}" "${MAX_LTR_PAIR_LENGTH_DIFF}" \
    >> "${TARGET}.sololtr_manifest.tsv"

echo "[OK] soloLTRseeker complete: ${TARGET} -> ${outdir}"
