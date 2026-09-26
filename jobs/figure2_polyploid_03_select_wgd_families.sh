#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# The input paths and commands are identical to the local runner.  jsub does
# not pass a script body through a wrapper reliably, so execute this file with:
#   jsub < jobs/figure2_polyploid_03_select_wgd_families.sh
exec bash jobs/figure2_polyploid_03_select_wgd_families_local.sh
