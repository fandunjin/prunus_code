#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Calculate corrected Figure 2 Ks from only newly generated JCVI/i-ADHoRe pairs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
JCVI_PAIRS="${JCVI_PAIRS:-${BASE}/01_jcvi_pairs/jcvi_syntenic_pairs.tsv}"
SELF_PAIRS="${SELF_PAIRS:-${BASE}/02_iadhore_self/self_syntenic_pairs.tsv}"
OUT="${OUT:-${BASE}/03_yn_kaks}"
THREADS="${THREADS:-12}"
MIN_AA="${MIN_AA:-50}"
MIN_CODONS="${MIN_CODONS:-50}"
MIN_4D_SITES="${MIN_4D_SITES:-4}"
MIN_REFERENCE_BLOCK_ANCHORS="${MIN_REFERENCE_BLOCK_ANCHORS:-20}"
MIN_CHROM_PAIR_ANCHORS="${MIN_CHROM_PAIR_ANCHORS:-20}"
MIN_CHROM_PAIR_FRACTION="${MIN_CHROM_PAIR_FRACTION:-0.005}"
KAKS_BIN="${KAKS_BIN:-${ROOT_DIR}/software/KaKs_Calculator3.0_linux/bin/KaKs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${BASE}/COMPLETE"
need "${JCVI_PAIRS}"
need "${SELF_PAIRS}"
need "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
for exe in mafft pal2nal.pl; do command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"; done
[[ -x "${KAKS_BIN}" ]] || die "KaKs_Calculator executable unavailable: ${KAKS_BIN}"

if [[ -e "${OUT}/COMPLETE" && "${FORCE}" != "1" ]]; then
  echo "[OK] Reusing corrected Figure 2 YN Ks results: ${OUT}"
  cat "${OUT}/ks_run_summary.tsv"
  exit 0
fi
if [[ -e "${OUT}" && "${FORCE}" == "1" ]]; then
  mv "${OUT}" "${OUT}.replaced.$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "${OUT}"

run_pairs() {
  local category="$1"
  local pair_table="$2"
  local category_out="${OUT}/${category}"
  echo "[INFO] ${category}: MAFFT -> PAL2NAL -> KaKs_Calculator 3.0 (YN)"
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py" \
    --analysis-dir "${ROOT_DIR}" --jcvi-dir "${BASE}" --pair-table "${pair_table}" \
    --outdir "${category_out}" --threads "${THREADS}" --pair-mode collinear_all_pairs \
    --min-aa "${MIN_AA}" --min-codons "${MIN_CODONS}" --min-4d-sites "${MIN_4D_SITES}" \
    --min-reference-block-anchors "${MIN_REFERENCE_BLOCK_ANCHORS}" \
    --min-chrom-pair-anchors "${MIN_CHROM_PAIR_ANCHORS}" \
    --min-chrom-pair-fraction "${MIN_CHROM_PAIR_FRACTION}" \
    --mafft-cmd mafft --pal2nal-cmd pal2nal.pl --kaks-cmd "${KAKS_BIN}" --kaks-method YN
  need "${category_out}/summary/all_merged_4dtv_kaks.tsv"
}

run_pairs jcvi "${JCVI_PAIRS}"
run_pairs self "${SELF_PAIRS}"

"${PYTHON_BIN}" - "${OUT}" <<'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])
rows = []
for category in ("jcvi", "self"):
    table = out / category / "summary" / "all_merged_4dtv_kaks.tsv"
    with table.open(encoding="utf-8", newline="") as handle:
        values = list(csv.DictReader(handle, delimiter="\t"))
    usable = 0
    for row in values:
        try:
            ks = float(row.get("Ks", ""))
            usable += int(0.0 < ks < 5.0)
        except ValueError:
            pass
    rows.append({"category": category, "aligned_pairs_with_kaks": len(values), "retained_0_lt_Ks_lt_5": usable, "source": table})
with (out / "ks_run_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY

{
  printf 'parameter\tvalue\n'
  printf 'input_contract\tcorrected_haploid_repartition_only_for_PZI_PSP_PPA\n'
  printf 'pair_sources\tJCVI_phase_specific_and_iADHoRe_real_cross_chromosome_anchorpoints\n'
  printf 'alignment\tMAFFT_protein_then_PAL2NAL_codon_backtranslation\n'
  printf 'ks_estimator\tKaKs_Calculator_3.0_YN\n'
  printf 'plot_filter\t0<Ks<5\n'
  printf 'min_aa\t%s\nmin_codons\t%s\nmin_4d_sites\t%s\n' "${MIN_AA}" "${MIN_CODONS}" "${MIN_4D_SITES}"
  printf 'min_reference_block_anchors\t%s\nmin_chrom_pair_anchors\t%s\nmin_chrom_pair_fraction\t%s\n' \
    "${MIN_REFERENCE_BLOCK_ANCHORS}" "${MIN_CHROM_PAIR_ANCHORS}" "${MIN_CHROM_PAIR_FRACTION}"
  printf 'threads\t%s\n' "${THREADS}"
} > "${OUT}/run_parameters.tsv"
printf 'Corrected Figure 2 YN Ks calculations completed.\n' > "${OUT}/COMPLETE"
echo "[OK] Corrected Figure 2 YN Ks: ${OUT}"
cat "${OUT}/ks_run_summary.tsv"
