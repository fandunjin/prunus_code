#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Build KaKs_Calculator 3.0 from the supplied source zip inside WSL.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
ZIP="${KAKS_ZIP:-${ROOT_DIR}/software/KaKs_Calculator3.0.zip}"
PREFIX="${KAKS_PREFIX:-${ROOT_DIR}/software/KaKs_Calculator3.0_linux}"
BUILD_ROOT="${BUILD_ROOT:-${ROOT_DIR}/software/.build_kaks_calculator3}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${ZIP}" ]] || die "Missing KaKs_Calculator source zip: ${ZIP}"
command -v unzip >/dev/null 2>&1 || die "unzip is unavailable on PATH"
command -v make >/dev/null 2>&1 || die "make is unavailable on PATH"
command -v g++ >/dev/null 2>&1 || die "g++ is unavailable on PATH"

rm -rf "${BUILD_ROOT}"
mkdir -p "${BUILD_ROOT}" "${PREFIX}/bin"
unzip -q "${ZIP}" -d "${BUILD_ROOT}"
SRC="${BUILD_ROOT}/KaKs_Calculator3.0/src"
[[ -d "${SRC}" ]] || die "Unexpected archive layout: ${SRC}"
make -C "${SRC}"
[[ -x "${SRC}/KaKs" ]] || die "KaKs build did not produce ${SRC}/KaKs"
cp "${SRC}/KaKs" "${PREFIX}/bin/KaKs"
chmod 755 "${PREFIX}/bin/KaKs"
printf '%s\n' "${PREFIX}/bin/KaKs" > "${PREFIX}/KAKS_EXECUTABLE"
echo "[OK] KaKs_Calculator 3.0 Linux binary: ${PREFIX}/bin/KaKs"
