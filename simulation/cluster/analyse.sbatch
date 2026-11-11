#!/bin/bash
# Collect and analyse the 2D Zeeman runs on the cluster: results/zeeman2d.dat, and summary.txt + figures in
# results/zeeman2d/ (plotting/zeeman2d_plots.py). Submitted automatically by submit_scan.sh with
# --dependency=afterany:<array job>, so it runs as soon as the whole scan has finished (also after failures:
# flagged runs are listed in summary.txt). By hand, from the repository root:
#     sbatch simulation/cluster/analyse.sbatch [--out DIR]       # DIR: results directory if not the default
#SBATCH --job-name=hcf2d_analyse
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=name.surname00@unibo.it

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source simulation/cluster/env.sh

OUT=""
if [ "${1:-}" = "--out" ]; then OUT="$2"; fi
echo "host $(hostname)  job ${SLURM_JOB_ID:-local}  start $(date)"
cd simulation/studies
if [ -n "$OUT" ]; then
    python zeeman2d.py collect --out "$OUT"
    python ../plotting/zeeman2d_plots.py --runs "$OUT"
else
    python zeeman2d.py collect
    python ../plotting/zeeman2d_plots.py
fi
echo "end $(date)"
