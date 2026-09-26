#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Validate the strict Figure 2a MUSCLE -> PAL2NAL -> CodeML interface on 10 pairs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
THREADS="${THREADS:-2}"

die() { echo "[ERROR] $*" >&2; exit 2; }

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
export PATH="${WGD_BIN}:${TOOLS_BIN}:/usr/bin:/bin"
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen

for exe in wgd muscle codeml; do
  command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
PREP="${BASE}/01_prepared_inputs"
TRACKS="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/01_catalog/tracks"
OUT="${BASE}/02_wgd_smoke_test"
FAMILIES_ALL="${PREP}/lower/PZI_A.families.tsv"
FAMILIES_TEST="${OUT}/PZI_A_GTR.first10.families.tsv"

[[ -s "${FAMILIES_ALL}" ]] || die "Missing input: ${FAMILIES_ALL}"
[[ -s "${TRACKS}/GTR.cds" ]] || die "Missing input: ${TRACKS}/GTR.cds"
[[ -s "${TRACKS}/PZI_A.cds" ]] || die "Missing input: ${TRACKS}/PZI_A.cds"

if [[ -e "${OUT}/COMPLETE" ]]; then
  die "Smoke test is already complete: ${OUT}"
fi
mkdir -p "${OUT}"
head -n 11 "${FAMILIES_ALL}" > "${FAMILIES_TEST}"
[[ "$(wc -l < "${FAMILIES_TEST}")" -eq 11 ]] || die "Failed to select 10 test families"

WGD_OUT="${OUT}/wgd_ksd"
WGD_TMP="${OUT}/tmp"
mkdir -p "${WGD_OUT}" "${WGD_TMP}"
KS_TSV="${WGD_OUT}/$(basename "${FAMILIES_TEST}").ks.tsv"

if [[ ! -s "${KS_TSV}" ]]; then
  wgd ksd \
    "${FAMILIES_TEST}" \
    "${TRACKS}/GTR.cds" \
    "${TRACKS}/PZI_A.cds" \
    --outdir "${WGD_OUT}" \
    --tmpdir "${WGD_TMP}" \
    --nthreads "${THREADS}" \
    --cds \
    --pairwise \
    --aligner muscle
else
  echo "[INFO] Reusing Ks table produced before the headless plotting failure: ${KS_TSV}"
fi

[[ -s "${KS_TSV}" ]] || die "Missing WGD Ks table: ${KS_TSV}"

"${WGD_BIN}/python3" - "${KS_TSV}" "${OUT}/smoke_test_summary.tsv" <<'PY'
import math
import sys
import pandas as pd

ks_path, summary_path = sys.argv[1:]
df = pd.read_csv(ks_path, sep="\t")
required = {"pair", "gene1", "gene2", "dS"}
missing = required - set(df.columns)
if missing:
    raise SystemExit("[ERROR] Missing WGD columns: {}".format(",".join(sorted(missing))))
values = pd.to_numeric(df["dS"], errors="coerce")
finite = values.map(math.isfinite)
positive = finite & (values > 0)
in_plot_range = positive & (values < 5)
summary = pd.DataFrame([{
    "families_input": 10,
    "rows_output": len(df),
    "finite_dS": int(finite.sum()),
    "positive_dS": int(positive.sum()),
    "dS_0_5": int(in_plot_range.sum()),
    "failed_or_nonfinite": int((~finite).sum()),
}])
summary.to_csv(summary_path, sep="\t", index=False)
if not finite.any():
    raise SystemExit("[ERROR] WGD produced no finite dS in the 10-family smoke test")
print(summary.to_string(index=False))
PY

{
  echo -e "component\tpath\tversion"
  echo -e "wgd\t$(command -v wgd)\t$(wgd --version 2>&1 | tail -n 1)"
  echo -e "muscle\t$(command -v muscle)\t$(muscle -version 2>&1 | head -n 1)"
  echo -e "codeml\t$(command -v codeml)\tPAML CodeML"
} > "${OUT}/tool_versions.tsv"

printf 'Strict Figure 2a WGD smoke test completed.\n' > "${OUT}/COMPLETE"
echo "[OK] Figure 2a WGD smoke test completed: ${OUT}"
echo "[OK] Ks table: ${KS_TSV}"
