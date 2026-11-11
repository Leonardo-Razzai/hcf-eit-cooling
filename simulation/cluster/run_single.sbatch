#!/bin/bash
# One simulation on the cluster (Slurm). Submit from the repository root, e.g.
#     sbatch simulation/cluster/run_single.sbatch zeeman2d.py check
#     sbatch simulation/cluster/run_single.sbatch zeeman2d.py run --scheme te01_bz --case B --s 1 --d1 0.4 --d2 0.4 --N 6
# Everything after the script name is a study script in simulation/studies and its arguments.
#
# Create the log folder once before the first submission: mkdir -p logs (Slurm needs it to exist).
# Adjust the #SBATCH lines to your needs; lines starting with "##SBATCH" are inactive.
# Times (4 cores, tmax = 1500/nu_r): N = 6 per axis ~1-2 h, N = 8 ~6-10 h, N = 10 ~1-2 days.
#SBATCH --job-name=hcf2d
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # BLAS threads (see env.sh)
#SBATCH --mem=4G                 # N <= 8: < 2 GB; N = 10: use 8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
##SBATCH --qos=debug             # max 15 min, no billing: use for check / smoke
##SBATCH --qos=long              # up to 72 h, low priority: N = 8, 10
##SBATCH --constraint=blade      # or matrix: ask your sector contact which to use
##SBATCH --account=oph           # only if your project requires it
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=name.surname00@unibo.it

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source simulation/cluster/env.sh

echo "host $(hostname)  job ${SLURM_JOB_ID:-local}  threads ${OPENBLAS_NUM_THREADS}  start $(date)"
SCRIPT="$1"; shift
cd simulation/studies
python "$SCRIPT" "$@"
echo "end $(date)"
