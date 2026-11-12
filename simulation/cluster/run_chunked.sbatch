#!/bin/bash
# A scan list in K chunks: array task k runs lines k, k+K, k+2K, ... one after the other (K = array size).
# Useful when the list is longer than the per-user job limits (normal QoS: 64 submitted, 8 running).
#     sbatch --array=1-8 simulation/cluster/run_chunked.sbatch simulation/cluster/scans/universality_fixedh.txt
# Finished runs are skipped (each run checks its output file), so resubmitting the same command resumes.
#SBATCH --job-name=hcf_chunk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
##SBATCH --constraint=matrix

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source simulation/cluster/env.sh

LIST="$1"; shift || true
K="${SLURM_ARRAY_TASK_COUNT:?run as an array job}"
k="${SLURM_ARRAY_TASK_ID}"
N=$(grep -cv '^\s*\(#\|$\)' "$LIST")
echo "host $(hostname)  job ${SLURM_ARRAY_JOB_ID}_${k} of ${K}  threads ${OPENBLAS_NUM_THREADS}  start $(date)"
for i in $(seq "$k" "$K" "$N"); do
    python simulation/cluster/run_list.py "$LIST" --line "$i" "$@" || echo "line $i failed"
done
echo "end $(date)"
