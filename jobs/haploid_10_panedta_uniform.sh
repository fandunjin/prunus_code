#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_ROOT="${TE_ROOT}/panEDTA_uniform"
WORK_ROOT="${OUT_ROOT}/work"
BACKUP_ROOT="${OUT_ROOT}/original_file_backup"
DIAG_ROOT="${OUT_ROOT}/suspect_family_diagnostics"
LOG_DIR="${ROOT_DIR}/logs"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
EDTA_ENV="${EDTA_ENV:-EDTA}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"

THREADS="${THREADS:-24}"
FL_COPY="${FL_COPY:-3}"
PHASE="${PHASE:-all}"                  # library or all
FORCE="${FORCE:-0}"                   # panEDTA -o: 0 resumes, 1 overwrites
RUN_DIAGNOSTICS="${RUN_DIAGNOSTICS:-1}"

HAPS=(hapA1 hapA2 hapB1 hapB2)

declare -A SUSPECT_TE=(
    [hapA1]="TE_00002654"
    [hapA2]="TE_00002302"
    [hapB1]="TE_00002357"
    [hapB2]="TE_00001842"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

require_command() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || die "Required command was not found: ${cmd}"
}

case "${PHASE}" in
    library|all) ;;
    *) die "PHASE must be 'library' or 'all', not '${PHASE}'" ;;
esac

case "${FORCE}" in
    0|1) ;;
    *) die "FORCE must be 0 or 1, not '${FORCE}'" ;;
esac

mkdir -p "${OUT_ROOT}" "${WORK_ROOT}" "${BACKUP_ROOT}" "${DIAG_ROOT}" "${LOG_DIR}"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${EDTA_CDS}" "CDS file used by panEDTA"

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${EDTA_ENV}"

require_command EDTA.pl
require_command rsync
require_command sha256sum

if command -v panEDTA.sh >/dev/null 2>&1; then
    PANEDTA_BIN="$(command -v panEDTA.sh)"
else
    PANEDTA_BIN="$(find "${CONDA_PREFIX}" -type f -name panEDTA.sh -print -quit 2>/dev/null || true)"
fi
[[ -n "${PANEDTA_BIN}" ]] || die "panEDTA.sh was not found in PATH or ${CONDA_PREFIX}"
require_file "${PANEDTA_BIN}" "panEDTA script"

GENOME_LIST="${OUT_ROOT}/pzi_four_haplotypes.list"
MANIFEST="${OUT_ROOT}/panEDTA_uniform_manifest.tsv"
CHECKSUMS="${BACKUP_ROOT}/original_files.sha256"
: > "${GENOME_LIST}"
: > "${CHECKSUMS}"

echo "=== Uniform panEDTA annotation for PZI haplotypes ==="
echo "EDTA root       : ${TE_ROOT}"
echo "Output root     : ${OUT_ROOT}"
echo "EDTA env        : ${EDTA_ENV}"
echo "panEDTA         : ${PANEDTA_BIN}"
echo "CDS             : ${EDTA_CDS}"
echo "Threads         : ${THREADS}"
echo "Full-copy cutoff: ${FL_COPY}"
echo "Phase           : ${PHASE}"
echo "Overwrite       : ${FORCE}"
echo

for hap in "${HAPS[@]}"; do
    src_dir="${TE_ROOT}/${hap}_edta"
    base="pzi_${hap}.fa"
    genome="${src_dir}/${base}"
    telib="${genome}.mod.EDTA.TElib.fa"
    gff="${genome}.mod.EDTA.TEanno.gff3"
    sum="${genome}.mod.EDTA.TEanno.sum"
    rmout="${src_dir}/${base}.mod.EDTA.anno/${base}.mod.EDTA.RM.out"

    require_file "${genome}" "${hap} genome"
    require_file "${telib}" "${hap} EDTA library"
    require_file "${gff}" "${hap} EDTA GFF3"
    require_file "${sum}" "${hap} EDTA summary"
    require_file "${rmout}" "${hap} EDTA RepeatMasker output"

    backup_dir="${BACKUP_ROOT}/${hap}"
    mkdir -p "${backup_dir}"
    for input in "${telib}" "${gff}" "${sum}" "${rmout}"; do
        output="${backup_dir}/$(basename "${input}")"
        if [[ ! -s "${output}" ]]; then
            cp -p "${input}" "${output}"
        fi
        sha256sum "${input}" >> "${CHECKSUMS}"
    done

    work_dir="${WORK_ROOT}/${hap}_edta"
    marker="${work_dir}/.source_edta_copy_complete"
    mkdir -p "${work_dir}"
    if [[ ! -e "${marker}" ]]; then
        echo "[COPY] ${src_dir} -> ${work_dir}"
        rsync -a "${src_dir}/" "${work_dir}/"
        touch "${marker}"
    else
        echo "[SKIP] Existing complete working copy: ${work_dir}"
    fi

    work_genome="${work_dir}/${base}"
    require_file "${work_genome}" "${hap} working genome"
    printf '%s\n' "${work_genome}" >> "${GENOME_LIST}"
done

echo
echo "Genome list: ${GENOME_LIST}"
cat "${GENOME_LIST}"

if [[ "${RUN_DIAGNOSTICS}" == "1" ]]; then
    if command -v seqkit >/dev/null 2>&1; then
        SUSPECT_FASTA="${DIAG_ROOT}/four_suspect_consensus.fa"
        : > "${SUSPECT_FASTA}"
        for hap in "${HAPS[@]}"; do
            te_id="${SUSPECT_TE[${hap}]}"
            lib="${TE_ROOT}/${hap}_edta/pzi_${hap}.fa.mod.EDTA.TElib.fa"
            out="${DIAG_ROOT}/${hap}.${te_id}.fa"
            seqkit grep -r -p "^${te_id}" "${lib}" > "${out}" || true
            if [[ -s "${out}" ]]; then
                if ! seqkit replace -p '^(.+)$' -r "${hap}|\$1" "${out}" >> "${SUSPECT_FASTA}"; then
                    echo "WARNING: failed to rename the ${hap} suspect-family sequence." >&2
                fi
            else
                echo "WARNING: suspect family ${te_id} was not found in ${lib}" >&2
            fi
        done

        if [[ -s "${SUSPECT_FASTA}" ]] && command -v minimap2 >/dev/null 2>&1; then
            if ! minimap2 -x asm20 -t "${THREADS}" "${SUSPECT_FASTA}" "${SUSPECT_FASTA}" \
                > "${DIAG_ROOT}/four_suspect_consensus.self.paf"; then
                echo "WARNING: minimap2 suspect-family self-alignment failed; continuing panEDTA." >&2
            fi
        fi
        if [[ -s "${SUSPECT_FASTA}" ]] && command -v TEsorter >/dev/null 2>&1; then
            if [[ ! -s "${SUSPECT_FASTA}.rexdb-plant.cls.tsv" || "${FORCE}" == "1" ]]; then
                if ! TEsorter "${SUSPECT_FASTA}" -db rexdb-plant -p "${THREADS}"; then
                    echo "WARNING: TEsorter diagnostics failed; continuing panEDTA." >&2
                fi
            fi
        fi
    else
        echo "WARNING: seqkit is unavailable; suspect-family diagnostics are skipped." >&2
    fi
fi

cd "${OUT_ROOT}"
PAN_ACTION=1
if [[ "${PHASE}" == "library" ]]; then
    PAN_ACTION=0
fi

echo
echo "=== Run panEDTA ==="
bash "${PANEDTA_BIN}" \
    -g "${GENOME_LIST}" \
    -c "${EDTA_CDS}" \
    -t "${THREADS}" \
    -f "${FL_COPY}" \
    -a "${PAN_ACTION}" \
    -o "${FORCE}"

EXPECTED_LIBRARY="${GENOME_LIST}.panEDTA.TElib.fa"
if [[ -s "${EXPECTED_LIBRARY}" ]]; then
    PAN_LIBRARY="${EXPECTED_LIBRARY}"
else
    PAN_LIBRARY="$(find "${OUT_ROOT}" -maxdepth 2 -type f -name '*.panEDTA.TElib.fa' -size +0c -print -quit)"
fi
[[ -n "${PAN_LIBRARY}" ]] || die "panEDTA completed but no non-empty panEDTA library was found"
require_file "${PAN_LIBRARY}" "Unified panEDTA library"

printf 'haplotype\tgenome\tpan_library\tpan_rmout\tteanno_gff3\tteanno_sum\tstatus\n' > "${MANIFEST}"
for hap in "${HAPS[@]}"; do
    work_dir="${WORK_ROOT}/${hap}_edta"
    base="pzi_${hap}.fa"
    genome="${work_dir}/${base}"
    pan_rmout="${genome}.mod.panEDTA.out"
    gff="${genome}.mod.EDTA.TEanno.gff3"
    sum="${genome}.mod.EDTA.TEanno.sum"

    status="library_only"
    if [[ "${PHASE}" == "all" ]]; then
        if [[ ! -s "${pan_rmout}" ]]; then
            pan_rmout="$(find "${work_dir}" -maxdepth 2 -type f -name '*.mod.panEDTA.out' -size +0c -print -quit)"
        fi
        [[ -n "${pan_rmout}" ]] || die "${hap} panEDTA RepeatMasker output was not found"
        require_file "${pan_rmout}" "${hap} panEDTA RepeatMasker output"
        require_file "${gff}" "${hap} unified TE annotation GFF3"
        require_file "${sum}" "${hap} unified TE annotation summary"
        status="complete"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${hap}" "${genome}" "${PAN_LIBRARY}" "${pan_rmout}" "${gff}" "${sum}" "${status}" \
        >> "${MANIFEST}"
done

echo
echo "=== panEDTA workflow completed ==="
echo "Unified library : ${PAN_LIBRARY}"
echo "Manifest        : ${MANIFEST}"
echo "Original backup : ${BACKUP_ROOT}"
echo "Working results : ${WORK_ROOT}"
echo "Diagnostics     : ${DIAG_ROOT}"
