#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
#
# Submit from the project checkout with:
#   jsub < jobs/psp_md5_check.sh

set -euo pipefail

ROOT="${ROOT:-${PROJECT_ROOT}/data/psp"
LIST="${PSP}/gsa_md5.list"

mkdir -p "${ROOT}/logs"

cat > "${LIST}" <<'EOF'
0edc46131c86a811179a8169abc15fc0  CRR2018654.bam
579cfb08633eab050992cc95f43992e2  CRR2018655_r1.fq.gz
0a4298e8c91ca00f5e5e4c2fc565a551  CRR2018655_r2.fq.gz
b6f880cc7ec147ad13225e3eb792e7b7  CRR2018656_r1.fq.gz
88b84e23893cb0db6f94607b4c0b239a  CRR2018656_r2.fq.gz
7e6fdd587c625e1dc4580cf077bb5704  CRR2018657_r1.fq.gz
045640d9e58324fbe26058df94451681  CRR2018657_r2.fq.gz
527ccf285b94f2086cc815e8dcbf3583  CRR2018658_r1.fq.gz
78866d6ef6c85a9883f60c6e0741ce8b  CRR2018658_r2.fq.gz
94b41b8acc7d490940583830c62d044d  CRR2018659_r1.fq.gz
34174e8e87ef6d55a92c14d38382f5b1  CRR2018659_r2.fq.gz
EOF

cd "${PSP}"
echo "=== md5 verification against GSA md5sum.txt ==="
md5sum -c "${LIST}" || true
echo
echo "Files whose line does not end in OK are incomplete or corrupted and must be re-downloaded."
