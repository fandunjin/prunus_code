#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
FINAL_PACKAGE_DIR="${ROOT_DIR}/final_haplotype_genomes"
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"

echo "=== Final PZI x=8 haplotype workflow ==="
echo
echo "1) Build evidence, score all chr01-32 four-chromosome candidates, and export source FASTA:"
echo "   jsub < jobs/haploid_01_build_phase_candidates.sh"
echo
echo "   Upstream source outputs (before C01-C08 renaming):"
echo "   ${FINAL_HAP_DIR}/FINAL_HAPLOTYPE_OVERVIEW.tsv"
echo "   ${FINAL_HAP_DIR}/FINAL_REVIEW_NOTES.tsv"
echo "   ${FINAL_HAP_DIR}/pzi_A1A2B1B2.rename_map.tsv"
echo "   ${FINAL_HAP_DIR}/pzi_hapA1.fa"
echo "   ${FINAL_HAP_DIR}/pzi_hapA2.fa"
echo "   ${FINAL_HAP_DIR}/pzi_hapB1.fa"
echo "   ${FINAL_HAP_DIR}/pzi_hapB2.fa"
echo
echo "2) Run Helixer for each final haplotype after FASTA extraction."
echo "   Recommended submission form:"
echo "   export HAP_SET=hapA1; jsub < jobs/haploid_03_helixer_final_hap.sh"
echo "   export HAP_SET=hapA2; jsub < jobs/haploid_03_helixer_final_hap.sh"
echo "   export HAP_SET=hapB1; jsub < jobs/haploid_03_helixer_final_hap.sh"
echo "   export HAP_SET=hapB2; jsub < jobs/haploid_03_helixer_final_hap.sh"
echo
echo "3) Postprocess each final haplotype annotation:"
echo "   export HAP_SET=hapA1; jsub < jobs/haploid_04_postprocess_final_hap.sh"
echo "   export HAP_SET=hapA2; jsub < jobs/haploid_04_postprocess_final_hap.sh"
echo "   export HAP_SET=hapB1; jsub < jobs/haploid_04_postprocess_final_hap.sh"
echo "   export HAP_SET=hapB2; jsub < jobs/haploid_04_postprocess_final_hap.sh"
echo
echo "4) Collect final downstream files after all four postprocess jobs finish:"
echo "   jsub < jobs/haploid_05_collect_final_hap_genomes.sh"
echo
echo "   Final minimal package:"
echo "   ${FINAL_PACKAGE_DIR}/"
echo
echo "5) Build the final C01-C08 chromosome-renamed downstream package:"
echo "   jsub < jobs/haploid_12c_rename_c_aliases.sh"
echo
echo "   Authoritative downstream package:"
echo "   ${HAP_RENAME_DIR}/"
echo "   Genome : ${HAP_RENAME_DIR}/genome/"
echo "   GFF    : ${HAP_RENAME_DIR}/gff/"
echo "   CDS    : ${HAP_RENAME_DIR}/cds/"
echo "   Protein: ${HAP_RENAME_DIR}/pep/"
echo
echo "Hard rule:"
echo "   Only x=8 is fixed. Quartet composition is selected by evidence-scored exact-cover logic."
