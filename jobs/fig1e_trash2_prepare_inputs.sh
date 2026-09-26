#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# fig1e_trash2_prepare_inputs.sh
# Run on the HPC login node: builds <FIG_ROOT>/inputs/genomes/<TARGET>.fa as
# symlinks (or decompression) from the existing HPC data locations.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/fig1e_satdna}"
INPUT_DIR="${INPUT_DIR:-${FIG_ROOT}/inputs/genomes}"

die() { echo "[ERROR] $*" >&2; exit 2; }

mkdir -p "${INPUT_DIR}"

link_input() {  # <source> <TARGET>
    local src="$1" target="$2"
    [[ -s "${src}" ]] || die "missing source: ${src}"
    ln -sf "${src}" "${INPUT_DIR}/${target}.fa"
    echo "  ${target}.fa -> ${src}"
}

decompress_input() {  # <source.gz> <TARGET>
    local src="$1" target="$2"
    [[ -s "${src}" ]] || die "missing source: ${src}"
    local dst="${INPUT_DIR}/${target}.fa"
    if [[ ! -s "${dst}" ]]; then
        gzip -dc "${src}" > "${dst}.tmp" && mv "${dst}.tmp" "${dst}"
    fi
    echo "  ${target}.fa (decompressed) -> ${src}"
}

HAP="${ROOT_DIR}/haploid_repartition/inputs"
DATA="${ROOT_DIR}/data"

echo "=== prepare pzi haplotypes ==="
link_input "${HAP}/pzi/genome/pzi_hapA1.fa" pziA1
link_input "${HAP}/pzi/genome/pzi_hapA2.fa" pziA2
link_input "${HAP}/pzi/genome/pzi_hapB1.fa" pziB1
link_input "${HAP}/pzi/genome/pzi_hapB2.fa" pziB2

echo "=== prepare psp haplotypes ==="
link_input "${HAP}/psp/genome/psp_hapA1.fa" pspA1
link_input "${HAP}/psp/genome/psp_hapA2.fa" pspA2
link_input "${HAP}/psp/genome/psp_hapB1.fa" pspB1
link_input "${HAP}/psp/genome/psp_hapB2.fa" pspB2

echo "=== prepare ppa subgenomes ==="
link_input "${HAP}/ppa/genome/ppa_hapA.fa" ppaA
link_input "${HAP}/ppa/genome/ppa_hapB.fa" ppaB

echo "=== prepare diploids and outgroups ==="
link_input "${DATA}/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna" ppe
link_input "${DATA}/gtr/GCA_018257905.1_V1_genomic.fna" gtr
decompress_input "${DATA}/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz" pcam
decompress_input "${DATA}/rch/GCF_002994745.2_RchiOBHm-V2_genomic.fna.gz" rch

echo "=== verify all 14 inputs ==="
ok=0
for t in pziA1 pziA2 pziB1 pziB2 pspA1 pspA2 pspB1 pspB2 ppaA ppaB ppe pcam gtr rch; do
    fa="${INPUT_DIR}/${t}.fa"
    if [[ -s "${fa}" ]]; then
        printf '  %-6s OK    %s bytes\n' "${t}" "$(stat -c %s "${fa}")"
        ok=$((ok + 1))
    else
        printf '  %-6s MISSING\n' "${t}"
    fi
done
[[ "${ok}" -eq 14 ]] || die "only ${ok}/14 inputs ready"
echo "=== inputs ready (${INPUT_DIR}) ==="
