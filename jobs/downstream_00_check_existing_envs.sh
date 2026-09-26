#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
REPORT_DIR="${ANALYSIS_DIR}/env_check"
REPORT="${REPORT_DIR}/software_env_check.tsv"
MISSING="${REPORT_DIR}/missing_software.tsv"
SUMMARY="${REPORT_DIR}/software_env_check_summary.txt"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

mkdir -p "${REPORT_DIR}" "${ANALYSIS_DIR}/logs"

if [ -s "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
else
    echo "WARNING: Conda init script not found: ${CONDA_SH}" >&2
fi

printf "environment\tsoftware\twhich\tversion_status\tversion\n" > "${REPORT}"

check_one() {
    local env_name="$1"
    local software="$2"
    local version_cmd="$3"
    local path
    local status
    local version
    local tmp_version
    local tmp_path
    local safe_name

    safe_name="${software//[^A-Za-z0-9_.-]/_}"
    tmp_path="${REPORT_DIR}/which.$$.$safe_name"
    which "${software}" > "${tmp_path}" 2>/dev/null || true
    path=""
    if [ -s "${tmp_path}" ]; then
        read -r path < "${tmp_path}"
    fi
    rm -f "${tmp_path}"
    if [ -z "${path}" ]; then
        printf "%s\t%s\t\tmissing\t\n" "${env_name}" "${software}" >> "${REPORT}"
        return
    fi
    tmp_version="${REPORT_DIR}/version.$$.$safe_name"
    set +e
    eval "${version_cmd}" > "${tmp_version}" 2>&1
    status=$?
    set -e
    version=""
    head -n 3 "${tmp_version}" > "${tmp_version}.head"
    while IFS= read -r line; do
        version="${version}${line};"
    done < "${tmp_version}.head"
    rm -f "${tmp_version}" "${tmp_version}.head"
    if [ "${status}" -eq 0 ]; then
        printf "%s\t%s\t%s\tok\t%s\n" "${env_name}" "${software}" "${path}" "${version}" >> "${REPORT}"
    else
        printf "%s\t%s\t%s\tversion_failed\t%s\n" "${env_name}" "${software}" "${path}" "${version}" >> "${REPORT}"
    fi
}

check_python_module() {
    local env_name="$1"
    local module="$2"
    local label="$3"
    local status
    local version
    local tmp_version
    local tmp_script
    local safe_label
    safe_label="${label//[^A-Za-z0-9_.-]/_}"
    tmp_version="${REPORT_DIR}/pyversion.$$.$safe_label"
    tmp_script="${REPORT_DIR}/pymodule.$$.$safe_label.py"
    cat > "${tmp_script}" <<EOF
import ${module}
print(getattr(${module}, "__version__", "import_ok"))
EOF
    set +e
    python "${tmp_script}" > "${tmp_version}" 2>&1
    status=$?
    set -e
    version=""
    head -n 3 "${tmp_version}" > "${tmp_version}.head"
    while IFS= read -r line; do
        version="${version}${line};"
    done < "${tmp_version}.head"
    rm -f "${tmp_version}" "${tmp_version}.head" "${tmp_script}"
    if [ "${status}" -eq 0 ]; then
        printf "%s\t%s\tpython-module\tok\t%s\n" "${env_name}" "${label}" "${version}" >> "${REPORT}"
    else
        printf "%s\t%s\t\tmissing\t%s\n" "${env_name}" "${label}" "${version}" >> "${REPORT}"
    fi
}

activate_env() {
    local env_name="$1"
    if command -v conda >/dev/null 2>&1; then
        set +e
        conda activate "${env_name}" >/dev/null 2>&1
        local status=$?
        set -e
        if [ "${status}" -ne 0 ]; then
            echo "WARNING: could not activate conda env '${env_name}', checking current PATH instead." >&2
        fi
    else
        echo "WARNING: conda command unavailable, checking current PATH for env '${env_name}'." >&2
    fi
}

echo "Checking base/current PATH"
check_one base bash "bash --version"
check_one base python "python --version"
check_one base python3 "python3 --version"
check_one base perl "perl -v"
check_one base java "java -version"
check_one base git "git --version"

activate_env "${PY_ENV:-haphic}"
check_one "${PY_ENV:-haphic}" python "python --version"
check_one "${PY_ENV:-haphic}" mafft "mafft --version"
check_one "${PY_ENV:-haphic}" trimal "trimal --version"
check_one "${PY_ENV:-haphic}" bmge "bmge -h"
check_one "${PY_ENV:-haphic}" pal2nal.pl "pal2nal.pl 2>&1"
check_one "${PY_ENV:-haphic}" ParaAT.pl "ParaAT.pl 2>&1"
check_one "${PY_ENV:-haphic}" KaKs_Calculator "KaKs_Calculator 2>&1"
check_one "${PY_ENV:-haphic}" seqmagick "seqmagick --version"
check_python_module "${PY_ENV:-haphic}" Bio biopython

activate_env "${BUSCO_ENV:-busco}"
check_one "${BUSCO_ENV:-busco}" busco "busco --version"

activate_env "${ORTHO_ENV:-orthofinder}"
check_one "${ORTHO_ENV:-orthofinder}" orthofinder "orthofinder --version"
check_one "${ORTHO_ENV:-orthofinder}" diamond "diamond version"
check_one "${ORTHO_ENV:-orthofinder}" mafft "mafft --version"
check_one "${ORTHO_ENV:-orthofinder}" iqtree3 "iqtree3 --version"
check_one "${ORTHO_ENV:-orthofinder}" iqtree2 "iqtree2 --version"
check_one "${ORTHO_ENV:-orthofinder}" iqtree "iqtree --version"

activate_env "${TREE_ENV:-phylo}"
check_one "${TREE_ENV:-phylo}" iqtree3 "iqtree3 --version"
check_one "${TREE_ENV:-phylo}" iqtree2 "iqtree2 --version"
check_one "${TREE_ENV:-phylo}" iqtree "iqtree --version"
check_one "${TREE_ENV:-phylo}" raxml-ng "raxml-ng --version"
check_one "${TREE_ENV:-phylo}" java "java -version"
check_one "${TREE_ENV:-phylo}" mcmctree "mcmctree 2>&1"
check_one "${TREE_ENV:-phylo}" baseml "baseml 2>&1"

activate_env "${CAFE_ENV:-cafe5}"
check_one "${CAFE_ENV:-cafe5}" cafe5 "cafe5 --version"

activate_env "${SYNTENY_ENV:-synteny}"
check_one "${SYNTENY_ENV:-synteny}" minimap2 "minimap2 --version"
check_one "${SYNTENY_ENV:-synteny}" nucmer "nucmer --version"
check_one "${SYNTENY_ENV:-synteny}" show-coords "show-coords -h"
check_one "${SYNTENY_ENV:-synteny}" python "python --version"
check_python_module "${SYNTENY_ENV:-synteny}" jcvi jcvi
check_one "${SYNTENY_ENV:-synteny}" MCScanX "MCScanX 2>&1"

activate_env "${REPEAT_ENV:-repeat}"
check_one "${REPEAT_ENV:-repeat}" EDTA.pl "EDTA.pl --help"
check_one "${REPEAT_ENV:-repeat}" LTR_retriever "LTR_retriever -h"
check_one "${REPEAT_ENV:-repeat}" RepeatModeler "RepeatModeler -h"
check_one "${REPEAT_ENV:-repeat}" RepeatMasker "RepeatMasker -h"
check_one "${REPEAT_ENV:-repeat}" gt "gt -version"

activate_env "${ANNOT_ENV:-annot}"
check_one "${ANNOT_ENV:-annot}" emapper.py "emapper.py --version"
check_one "${ANNOT_ENV:-annot}" interproscan.sh "interproscan.sh -version"

activate_env "${NETWORK_ENV:-network}"
check_one "${NETWORK_ENV:-network}" java "java -version"
check_one "${NETWORK_ENV:-network}" python "python --version"
check_one "${NETWORK_ENV:-network}" run_hyde.py "run_hyde.py -h"
check_one "${NETWORK_ENV:-network}" phyde "phyde -h"
check_one "${NETWORK_ENV:-network}" Rscript "Rscript --version"

awk -F '\t' 'NR == 1 || $4 == "missing" || $4 == "version_failed"' "${REPORT}" > "${MISSING}"
{
    echo "Software environment check report: ${REPORT}"
    echo "Missing/version-failed report: ${MISSING}"
    echo
    echo "Missing or version-failed entries:"
    tail -n +2 "${MISSING}" || true
} > "${SUMMARY}"

cat "${SUMMARY}"
