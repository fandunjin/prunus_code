#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Create or repair the isolated strict Figure 2a analysis environment.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
ENV_NAME="${ENV_NAME:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"

die() { echo "[ERROR] $*" >&2; exit 2; }
command -v conda >/dev/null 2>&1 || die "conda is unavailable on PATH"

CONDA_BASE="$(conda info --base)"
[[ -s "${CONDA_BASE}/etc/profile.d/conda.sh" ]] || die "Cannot find conda.sh under ${CONDA_BASE}"
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

# The official WGD suite is installed from the heche-psb/wgd v2 source below.
# The similarly named Conda package has a different CLI and must not be used.
# Python dependencies are installed from the upstream requirements.txt. The
# pre-existing tools environment supplies only non-Python executables; keeping
# its libraries out of this environment avoids Conda/pip ABI conflicts.
PACKAGES=(
    python=3.8 pip setuptools wheel
)
# Ignore stale global mirrors for this installation only. Override with
# CONDA_CHANNEL_ARGS if the cluster requires an institutional mirror.
if [[ -n "${CONDA_CHANNEL_ARGS:-}" ]]; then
    read -r -a CHANNEL_ARGS <<< "${CONDA_CHANNEL_ARGS}"
else
    CHANNEL_ARGS=(
        -c https://conda.anaconda.org/bioconda
        -c https://conda.anaconda.org/conda-forge
    )
fi
if conda env list | awk 'NR > 2 {print $1}' | grep -Fxq "${ENV_NAME}"; then
    conda install -y -n "${ENV_NAME}" --override-channels "${CHANNEL_ARGS[@]}" "${PACKAGES[@]}"
else
    conda create -y -n "${ENV_NAME}" --override-channels "${CHANNEL_ARGS[@]}" "${PACKAGES[@]}"
fi

conda env list | awk 'NR > 2 {print $1}' | grep -Fxq "${TOOLS_ENV}" || die "Missing tools environment: ${TOOLS_ENV}"
ENV_PREFIX="${CONDA_BASE}/envs/${ENV_NAME}"
TOOLS_PREFIX="${CONDA_BASE}/envs/${TOOLS_ENV}"
ENV_PYTHON="${ENV_PREFIX}/bin/python"
TOOLS_BIN="${TOOLS_PREFIX}/bin"
[[ -x "${ENV_PYTHON}" ]] || die "Cannot resolve Python from clean environment: ${ENV_PREFIX}"
[[ -x "${TOOLS_BIN}/blastp" ]] || die "Cannot resolve blastp from tools environment: ${TOOLS_ENV}"
export WGD2_TOOLS_BIN="${TOOLS_BIN}"
export PATH="${ENV_PREFIX}/bin:${WGD2_TOOLS_BIN}:${PATH}"

# The upstream importlib-metadata pin predates the modern setuptools API.
"${ENV_PYTHON}" -m pip install --no-cache-dir "setuptools<58"

WGD_ROOT="${ROOT_DIR}/software/wgd_v2"
[[ -f "${WGD_ROOT}/setup.py" ]] || die "Missing official wgd_v2 source: ${WGD_ROOT}"
[[ -f "${WGD_ROOT}/requirements.txt" ]] || die "Missing official wgd_v2 requirements: ${WGD_ROOT}"
WGD_COMMIT="$(git -C "${WGD_ROOT}" rev-parse HEAD 2>/dev/null || true)"
[[ -n "${WGD_COMMIT}" ]] || die "Cannot record the wgd_v2 source commit"
EXPECTED_WGD_COMMIT="${EXPECTED_WGD_COMMIT:-d2c06a543bc272b2cb5dc469f600380e7dc0771a}"
[[ "${WGD_COMMIT}" == "${EXPECTED_WGD_COMMIT}" ]] || die "Unexpected wgd_v2 commit ${WGD_COMMIT}; expected ${EXPECTED_WGD_COMMIT}"

# Upstream's documented installation order. The pinned stack is required
# because wgd v2 resolves its installed requirements through pkg_resources.
"${ENV_PYTHON}" -m pip install --no-cache-dir "numpy==1.19.0"
"${ENV_PYTHON}" -m pip install --no-cache-dir -r "${WGD_ROOT}/requirements.txt"
"${ENV_PYTHON}" -m pip install --no-deps --no-build-isolation --upgrade "${WGD_ROOT}"

WGD_VERSION="$("${ENV_PYTHON}" -c 'import importlib.metadata as m; print(m.version("wgd"))' 2>/dev/null || true)"
[[ "${WGD_VERSION}" == "2.0.38" ]] || die "Expected official wgd 2.0.38, found ${WGD_VERSION:-missing}"
[[ "$("${ENV_PYTHON}" -c 'import Bio; print(Bio.__version__)')" == "1.76" ]] || die "Expected Biopython 1.76"
[[ "$("${ENV_PYTHON}" -c 'import numpy; print(numpy.__version__)')" == "1.19.0" ]] || die "Expected NumPy 1.19.0"
WGD_ENTRY="${ENV_PREFIX}/bin/wgd"
[[ -n "${WGD_ENTRY}" ]] || die "Official wgd executable is unavailable after installation"
[[ -x "${WGD_ENTRY}" ]] || die "Official wgd executable is unavailable after installation"
"${ENV_PYTHON}" -c 'import cli; print(cli.__file__)' >/dev/null 2>&1 || die "wgd entry module cannot be imported"
"${WGD_ENTRY}" --help >/dev/null 2>&1 || die "wgd --help failed"
"${WGD_ENTRY}" dmd --help >/dev/null 2>&1 || die "wgd dmd --help failed"
"${WGD_ENTRY}" ksd --help >/dev/null 2>&1 || die "wgd ksd --help failed"

MANIFEST="${ROOT_DIR}/software/panel_a_reference_env.tsv"
printf 'tool\tpath\tversion\n' > "${MANIFEST}"
for tool in python wgd blastp makeblastdb mcl muscle pal2nal.pl codeml; do
    path="$(command -v "${tool}" 2>/dev/null || true)"
    version="unknown"
    case "${tool}" in
        python) version="$("${ENV_PYTHON}" --version 2>&1 | head -n 1)" ;;
        wgd) version="${WGD_VERSION}" ;;
        blastp|makeblastdb) version="$("${tool}" -version 2>&1 | head -n 1 || true)" ;;
        mcl) version="$(mcl --version 2>&1 | head -n 1 || true)" ;;
        muscle) version="$(muscle -version 2>&1 | head -n 1 || true)" ;;
        pal2nal.pl) version="$(pal2nal.pl 2>&1 | head -n 1 || true)" ;;
        codeml) version="PAML CodeML (version will be captured from formal run output)" ;;
    esac
    printf '%s\t%s\t%s\n' "${tool}" "${path:-NA}" "${version:-unknown}" >> "${MANIFEST}"
done

grep -q $'^muscle\t.*3\.8\.1551' "${MANIFEST}" || die "Installed MUSCLE is not 3.8.1551; see ${MANIFEST}"
printf 'wgd_source\t%s\n' "${WGD_ROOT}" >> "${MANIFEST}"
printf 'wgd_version\t%s\n' "${WGD_VERSION}" >> "${MANIFEST}"
printf 'wgd_entry\t%s\n' "${WGD_ENTRY}" >> "${MANIFEST}"
printf 'wgd_commit\t%s\n' "${WGD_COMMIT}" >> "${MANIFEST}"
printf 'tools_environment\t%s\n' "${TOOLS_ENV}" >> "${MANIFEST}"
printf 'tools_bin\t%s\n' "${TOOLS_BIN}" >> "${MANIFEST}"
echo "[OK] Strict Figure 2a environment ready: ${ENV_NAME}"
echo "[OK] Environment manifest: ${MANIFEST}"
