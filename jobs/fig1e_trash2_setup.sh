#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# fig1e_trash2_setup.sh
# Run ONCE on the HPC login node (bash jobs/fig1e_trash2_setup.sh).
# Creates the fig1e conda env (R + TRASH_2 dependencies), installs TRASH_2
# (patched bundle preferred) and verifies the environment.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/fig1e_satdna}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BASE_PY="${BASE_PY:-${CONDA_ROOT}/bin/python3}"

die() { echo "[ERROR] $*" >&2; exit 2; }

[[ -s "${CONDA_SH}" ]] || die "missing ${CONDA_SH}"
source "${CONDA_SH}"

if conda env list | awk '{print $1}' | grep -qx "fig1e"; then
    echo "[SKIP] conda env fig1e exists"
else
    echo "=== create conda env fig1e ==="
    conda create -y -q -n fig1e -c conda-forge -c bioconda \
        r-base mafft hmmer gcc_linux-64 gxx_linux-64 make zlib libcurl libxml2 \
        r-rcpp r-stringr r-stringdist r-seqinr r-doparallel r-getopt r-ape \
        r-biocmanager r-igraph r-domc bioconductor-biostrings
fi

conda activate fig1e

echo "=== install msa (Bioconductor, compiles once) ==="
Rscript -e 'if (!requireNamespace("msa", quietly=TRUE)) BiocManager::install("msa", update=FALSE, ask=FALSE)'

echo "=== install BCT (CRAN, used for CTW compression tracks) ==="
Rscript -e 'if (!requireNamespace("BCT", quietly=TRUE)) install.packages("BCT", repos="https://cloud.r-project.org")'

echo "=== install TRASH_2 ==="
mkdir -p "${FIG_ROOT}/tools"
if [[ -d "${FIG_ROOT}/tools/TRASH_2/src" ]]; then
    echo "[SKIP] ${FIG_ROOT}/tools/TRASH_2 already present"
elif [[ -s "${FIG_ROOT}/tools/TRASH_2_patched.tar.gz" ]]; then
    tar xzf "${FIG_ROOT}/tools/TRASH_2_patched.tar.gz" -C "${FIG_ROOT}/tools"
    echo "[OK] extracted patched bundle"
else
    git clone https://github.com/vlothec/TRASH_2.git "${FIG_ROOT}/tools/TRASH_2"
fi

TRASH_MAIN="${FIG_ROOT}/tools/TRASH_2/src/main.R"
[[ -s "${TRASH_MAIN}" ]] || die "missing ${TRASH_MAIN}"

if ! grep -q "doMC" "${TRASH_MAIN}"; then
    echo "=== applying doMC patch ==="
    [[ -x "${BASE_PY}" ]] || die "missing python3 for patching: ${BASE_PY}"
    "${BASE_PY}" - "${TRASH_MAIN}" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
s = p.read_text()
old = """  cl <- makeCluster(cmd_arguments$cores_no,
                    # outfile="",
                    homogeneous = TRUE)
  clusterEvalQ(cl, .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd()))))
  clusterEvalQ(cl, sink())
  registerDoParallel(cl)"""
new = """  if ((Sys.info()["sysname"] == "Linux") && requireNamespace("doMC", quietly = TRUE)) {
    doMC::registerDoMC(cmd_arguments$cores_no)
    cl <- NULL
  } else {
    cl <- makeCluster(cmd_arguments$cores_no,
                      # outfile="",
                      homogeneous = TRUE)
    clusterEvalQ(cl, .libPaths(c(.libPaths(), gsub("src", "R_libs", getwd()))))
    clusterEvalQ(cl, sink())
    registerDoParallel(cl)
  }"""
assert old in s, "cluster block not found"
s = s.replace(old, new, 1)
s = s.replace("  stopCluster(cl)", "  if (!is.null(cl)) stopCluster(cl)", 1)
p.write_text(s)
print("patched", p)
PY
fi
grep -q "doMC" "${TRASH_MAIN}" || die "doMC patch not present"

echo "=== verify R packages ==="
Rscript -e 'pkgs <- c("stringr","stringdist","seqinr","doParallel","getopt","ape","Biostrings","msa","doMC","BCT"); ok <- sapply(pkgs, requireNamespace, quietly=TRUE); print(ok); quit(status = !all(ok))'

command -v mafft >/dev/null || die "mafft missing"
command -v nhmmer >/dev/null || die "nhmmer missing"
command -v Rscript >/dev/null || die "Rscript missing"

echo "=== setup done; TRASH_2: ${FIG_ROOT}/tools/TRASH_2 ==="
