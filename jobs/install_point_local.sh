#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Build the official POInT source tree already cloned under software/POInT.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SOURCE_DIR="${SOURCE_DIR:-${ROOT_DIR}/software/POInT}"
THREADS="${THREADS:-12}"
POINT_OPENMP="${POINT_OPENMP:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -f "${SOURCE_DIR}/configure.pl" ]] || die "Missing official POInT source: ${SOURCE_DIR}"
for tool in perl make gcc g++; do
    command -v "${tool}" >/dev/null 2>&1 || die "Required build tool unavailable: ${tool}"
done

(
    cd "${SOURCE_DIR}"
    make clean >/dev/null 2>&1 || true
    if [[ "${POINT_OPENMP}" == "1" ]]; then
        perl configure.pl -omp
    else
        perl configure.pl
    fi
    make -j "${THREADS}"
)

MANIFEST="${SOURCE_DIR}/BUILD_MANIFEST.tsv"
printf 'field\tvalue\n' > "${MANIFEST}"
printf 'commit\t%s\n' "$(git -C "${SOURCE_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)" >> "${MANIFEST}"
printf 'threads\t%s\n' "${THREADS}" >> "${MANIFEST}"
if [[ "${POINT_OPENMP}" == "1" ]]; then
    printf 'configure\tperl configure.pl -omp\n' >> "${MANIFEST}"
else
    printf 'configure\tperl configure.pl\n' >> "${MANIFEST}"
fi
printf 'openmp\t%s\n' "${POINT_OPENMP}" >> "${MANIFEST}"
printf 'local_source_fixes\tvirtual Branch/Tree destructors; delete[] for Phylo_Matrix matrix_params rows\n' >> "${MANIFEST}"
printf 'executables\t%s\n' "$(find "${SOURCE_DIR}" -maxdepth 1 -type f -perm -111 -printf '%f,' | sed 's/,$//')" >> "${MANIFEST}"
[[ -x "${SOURCE_DIR}/POInT" ]] || die "POInT build did not produce ${SOURCE_DIR}/POInT"
[[ -x "${SOURCE_DIR}/POInT_ances_order" ]] || die "POInT build did not produce ${SOURCE_DIR}/POInT_ances_order"
touch "${SOURCE_DIR}/BUILD_COMPLETE"
echo "[OK] Official POInT built: ${SOURCE_DIR}"
echo "[OK] Build manifest: ${MANIFEST}"
