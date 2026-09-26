#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Strict reference-method Figure 2a software and input preflight.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
CATALOG="${CATALOG:-${BASE}/04_fractbias_8track_reference_exact/01_catalog/tracks}"
OUT="${OUT:-${BASE}/06_panel_a_reference_exact/00_preflight}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VALIDATOR="${ROOT_DIR}/scripts/fig2_panel_a_validate_inputs.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

[[ ! -e "${OUT}/COMPLETE" ]] || die "Preflight is already complete: ${OUT}"
require_file "${VALIDATOR}"
mkdir -p "${OUT}"

"${PYTHON_BIN}" "${VALIDATOR}" --tracks "${CATALOG}" --out "${OUT}/input_manifest.tsv"

TOOL_TABLE="${OUT}/tool_preflight.tsv"
printf 'tool\tpath\tversion\tstatus\tnote\n' > "${TOOL_TABLE}"
missing=0
version_fail=0
interface_fail=0

resolve_tool() {
    local name="$1"
    local path=""
    if command -v "${name}" >/dev/null 2>&1; then
        path="$(command -v "${name}")"
    elif [[ "${name}" == "i-adhore" && -x "${ROOT_DIR}/software/i-ADHoRe/install/bin/i-adhore" ]]; then
        path="${ROOT_DIR}/software/i-ADHoRe/install/bin/i-adhore"
    fi
    printf '%s' "${path}"
}

record_tool() {
    local name="$1"
    local path version status note
    path="$(resolve_tool "${name}")"
    if [[ -z "${path}" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\n' "${name}" NA NA MISSING required >> "${TOOL_TABLE}"
        missing=$((missing + 1))
        return
    fi
    version="unknown"
    status="PASS"
    note="required"
    case "${name}" in
        python3) version="$(${path} --version 2>&1 | head -n 1)" ;;
        wgd)
            version="$(python3 -c 'import importlib.metadata as m; print(m.version("wgd"))' 2>/dev/null || true)"
            if [[ "${version}" != "2.0.38" ]]; then
                status="WRONG_VERSION"
                note="reference requires official wgd v2 2.0.38"
                version_fail=$((version_fail + 1))
            elif ! "${path}" dmd --help >/dev/null 2>&1 || ! "${path}" ksd --help >/dev/null 2>&1; then
                status="INTERFACE_FAIL"
                note="official wgd dmd/ksd interfaces unavailable"
                interface_fail=$((interface_fail + 1))
            fi
            ;;
        blastp|makeblastdb) version="$(${path} -version 2>&1 | head -n 1 || true)" ;;
        mcl) version="$(${path} --version 2>&1 | head -n 1 || true)" ;;
        muscle) version="$(${path} -version 2>&1 | head -n 1 || true)" ;;
        pal2nal.pl) version="$(${path} 2>&1 | head -n 1 || true)" ;;
        codeml) version="PAML CodeML; version recorded from formal run output" ;;
        i-adhore) version="$(${path} --version 2>&1 | head -n 1 || true)" ;;
    esac
    [[ -n "${version}" ]] || version="unknown"
    if [[ "${name}" == "muscle" && "${version}" != *"3.8.1551"* ]]; then
        status="WRONG_VERSION"
        note="reference requires MUSCLE 3.8.1551"
        version_fail=$((version_fail + 1))
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "${name}" "${path}" "${version}" "${status}" "${note}" >> "${TOOL_TABLE}"
}

for tool in python3 wgd blastp makeblastdb mcl muscle pal2nal.pl codeml i-adhore; do
    record_tool "${tool}"
done

cat > "${OUT}/method_contract.txt" <<'EOF'
Panel a reference-method contract:
- WGD/paranome and Ks interface: official wgd v2.0.38 (heche-psb/wgd),
  recorded by source commit and environment manifest;
- homolog search: all-against-all BLASTP, E <= 1e-10;
- family clustering: MCL (inflation and any unreported defaults must be recorded);
- A-B collinearity: i-ADHoRe 3.0, never JCVI anchors;
- protein alignment: MUSCLE 3.8.1551;
- codon ML: PAML CodeML, Nielsen-Yang maximum-likelihood model;
- CDS records whose length is not a multiple of three are excluded before
  codon alignment and are counted in the final filtering table;
- upper plot: PZI, PSP and PPA A-B collinear pairs;
- lower plot: GTR-to-PZI/PSP/PPA/PCAM/PPE ortholog context;
- Ks is a relative divergence signal, not a direct date or proof of a shared event.
EOF

if (( missing > 0 || version_fail > 0 || interface_fail > 0 )); then
    printf 'status\tBLOCKED\nmissing_tools\t%d\nversion_failures\t%d\ninterface_failures\t%d\n' "${missing}" "${version_fail}" "${interface_fail}" > "${OUT}/preflight_status.tsv"
    echo "[ERROR] Panel a preflight blocked: missing=${missing}; wrong_version=${version_fail}; interface_fail=${interface_fail}" >&2
    echo "[INFO] See ${TOOL_TABLE}" >&2
    exit 2
fi

printf 'status\tPASS\nmissing_tools\t0\nversion_failures\t0\ninterface_failures\t0\n' > "${OUT}/preflight_status.tsv"
touch "${OUT}/COMPLETE"
echo "[OK] Strict Figure 2a preflight completed: ${OUT}"
