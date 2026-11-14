# Release validation

Checked against the manuscript analysis outputs on 2026-11-14.

| Check | Result |
| --- | --- |
| Representative-locus family analysis | Summary and eight result tables match; 53 retained families |
| POInT model comparison | All five model rows match, including likelihoods, parameter counts, AIC and BIC |
| Regional retention profiles | All 3,822 window rows match |
| GO meta-analysis | Five input/output tables match |
| CAFE input preparation | Full and diploid-control inputs each match all 20,239 rows |
| Satellite representatives | Five sequence/metadata files match; 168 alignment columns |
| Genome-size decomposition | All 22 component rows match |
| Pairwise chromosome inputs | Three BED and three protein files match byte-for-byte |
| Figure 2 source blocks | All 341 rows match, including the 29 non-adjacent PZI–PPA blocks; adjacent pairs total 312 |
| JCVI display coordinates | All 80 chromosome-layout rows match |
| Automated tests | Six tests pass: locus identity/selection, multiple testing, partial conjunction and repeat-overlap priorities |
| Command-line interfaces | All 46 executable scripts respond to `--help`; three shared library modules are importable |
| Static checks | All Python files parse; no undefined-name/unused-import diagnostics; no whitespace errors |
| Repository cleanup | 49 current analysis/helper scripts, down from 673 historical scripts; no duplicate version files |
| Privacy review | No personal absolute paths, account/host identifiers, scheduler directives or credentials in the current files |

The release comparisons rerun the affected data-preparation and summary calculations.
Genome assembly, full OrthoFinder inference and all codeml fits were not rerun as
part of this code cleanup. Their analysis commands and parameters are documented
in [commands.md](commands.md).
