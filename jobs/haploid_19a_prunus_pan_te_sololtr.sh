#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/prunus_pan_te}"
FOLLOWUP_ROOT="${FOLLOWUP_ROOT:-${OUT_ROOT}/followup}"
TARGET_CONFIG="${TARGET_CONFIG:-${ROOT_DIR}/scripts/prunus_pan_te_targets.sh}"
SOLOLTR_DIR="${SOLOLTR_DIR:-${ROOT_DIR}/software/soloLTRseeker/modules}"
SOLOLTR_BIN="${SOLOLTR_BIN:-}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SOLOLTR_ENV="${SOLOLTR_ENV:-soloLTR}"
THREADS="${THREADS:-8}"
MIN_LTR_LENGTH="${MIN_LTR_LENGTH:-100}"
MAX_LTR_PAIR_LENGTH_DIFF="${MAX_LTR_PAIR_LENGTH_DIFF:-50}"

[ -n "${TARGET}" ] || { echo "ERROR: set TARGET" >&2; exit 1; }
[ -s "${TARGET_CONFIG}" ] || { echo "ERROR: missing ${TARGET_CONFIG}" >&2; exit 1; }
[ -s "${CONDA_SH}" ] || { echo "ERROR: missing ${CONDA_SH}" >&2; exit 1; }
source "${TARGET_CONFIG}"
source "${CONDA_SH}"
conda activate "${SOLOLTR_ENV}"

# Resolve soloLTRseeker after environment activation.  The upstream repository
# has used both modules/soloLTRseeker and top-level soloLTRseeker layouts.
if [ -z "${SOLOLTR_BIN}" ]; then
    for candidate in \
        "${SOLOLTR_DIR}/soloLTRseeker" \
        "${ROOT_DIR}/software/soloLTRseeker/soloLTRseeker" \
        "${ROOT_DIR}/software/soloLTRseeker/bin/soloLTRseeker"; do
        if [ -s "${candidate}" ]; then
            SOLOLTR_BIN="${candidate}"
            break
        fi
    done
fi
if [ -z "${SOLOLTR_BIN}" ]; then
    SOLOLTR_BIN="$(command -v soloLTRseeker || true)"
fi
if [ -z "${SOLOLTR_BIN}" ] || [ ! -s "${SOLOLTR_BIN}" ]; then
    cat >&2 <<EOF
ERROR: soloLTRseeker was not found.
Install it at ${ROOT_DIR}/software/soloLTRseeker or set SOLOLTR_BIN to its executable path.
Expected upstream layout: ${ROOT_DIR}/software/soloLTRseeker/modules/soloLTRseeker
EOF
    exit 1
fi
if [ ! -x "${SOLOLTR_BIN}" ]; then
    echo "ERROR: soloLTRseeker is not executable: ${SOLOLTR_BIN}" >&2
    echo "Run: chmod u+x ${SOLOLTR_BIN}" >&2
    exit 1
fi
set_pan_te_target "${TARGET}"
for exe in bedtools blastn cd-hit-est water TEsorter; do
    command -v "${exe}" >/dev/null 2>&1 || { echo "ERROR: ${exe} was not found" >&2; exit 1; }
done

genome="${OUT_ROOT}/discovery/${TARGET}/${TARGET}.fa"
intact_gff="${genome}.mod.EDTA.intact.gff3"
outdir="${FOLLOWUP_ROOT}/sololtr/${TARGET}"
for path in "${genome}" "${intact_gff}"; do
    [ -s "${path}" ] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
mkdir -p "${outdir}"
cd "${outdir}"

# soloLTRseeker parses the complete FASTA header and its output-prefix logic
# expects a simple *.fa basename.  Build a compute-node-local FASTA that keeps
# only the sequence identifier before the first whitespace.
sololtr_genome="${outdir}/${TARGET}.fa"
if [ ! -s "${sololtr_genome}" ]; then
    tmp_genome="${sololtr_genome}.tmp.$$"
    awk '/^>/{sub(/[[:space:]].*$/, "")}1' "${genome}" > "${tmp_genome}"
    [ -s "${tmp_genome}" ] || { echo "ERROR: failed to create ${tmp_genome}" >&2; exit 1; }
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
        "${intact_gff}" "${sololtr_genome}"
fi

solo_file="$(find . -type f -name '*_soloLTR.gff3' -size +0c -printf '%T@\t%p\n' \
    | sort -nr | head -n 1 | cut -f2-)"
length_file="$(find . -type f -name '*_LTRRT_lengths.txt' -size +0c -printf '%T@\t%p\n' \
    | sort -nr | head -n 1 | cut -f2-)"
[ -n "${solo_file}" ] || { echo "ERROR: no non-empty soloLTR GFF found recursively under ${outdir}" >&2; exit 1; }
[ -n "${length_file}" ] || { echo "ERROR: no non-empty LTRRT length table found recursively under ${outdir}" >&2; exit 1; }
ln -sfn "$(readlink -f "${solo_file}")" "${TARGET}.soloLTR.gff3"
ln -sfn "$(readlink -f "${length_file}")" "${TARGET}.LTRRT_lengths.txt"

printf "target\tgenome\tintact_gff\tsolo_gff\tltrrt_lengths\tmin_ltr_length\tmax_pair_length_difference\n" > "${TARGET}.sololtr_manifest.tsv"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${TARGET}" "${genome}" "${intact_gff}" "$(readlink -f "${TARGET}.soloLTR.gff3")" \
    "$(readlink -f "${TARGET}.LTRRT_lengths.txt")" "${MIN_LTR_LENGTH}" "${MAX_LTR_PAIR_LENGTH_DIFF}" \
    >> "${TARGET}.sololtr_manifest.tsv"

echo "=== soloLTRseeker completed: ${TARGET} ==="
