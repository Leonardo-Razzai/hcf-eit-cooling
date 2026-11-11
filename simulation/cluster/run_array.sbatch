#!/bin/bash
# Parameter scan as a Slurm job array: one array task per line of a scan list
# (written by "python zeeman2d.py list ..." into simulation/cluster/scans/). Submit from the repository root:
#
#     LIST=simulation/cluster/scans/zeeman2d_main.txt
#     N=$(grep -cv '^\s*\(#\|$\)' $LIST)
#     sbatch --array=1-$N%20 simulation/cluster/run_array.sbatch $LIST
#
# "%20" limits the number of tasks running at the same time (be fair to other users).
# Finished runs (their JSON file exists) are skipped, so the same command can be resubmitted
# after a failure or a timeout.
#SBATCH --job-name=hcf2d_scan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
##SBATCH --qos=long
##SBATCH --constraint=blade
##SBATCH --account=oph
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=name.surname00@unibo.it

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source simulation/cluster/env.sh

LIST="$1"; shift || true
echo "host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-?}  threads ${OPENBLAS_NUM_THREADS}  start $(date)"
python simulation/cluster/run_list.py "$LIST" --line "${SLURM_ARRAY_TASK_ID:?run as an array job}" "$@"
echo "end $(date)"
