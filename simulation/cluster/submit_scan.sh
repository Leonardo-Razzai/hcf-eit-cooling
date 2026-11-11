#!/bin/bash
# Submit a scan list as a job array and, after it, the analysis job (collect + summary + figures).
# Run on the frontend from the repository root:
#     bash simulation/cluster/submit_scan.sh simulation/cluster/scans/zeeman2d_main.txt [sbatch options for the array]
# e.g. for the convergence list:
#     bash simulation/cluster/submit_scan.sh simulation/cluster/scans/zeeman2d_convergence.txt --qos=long --time=72:00:00 --mem=8G --cpus-per-task=8
# The job ids are appended to logs/submitted.txt (read by remote.sh status).
set -euo pipefail
LIST="$1"; shift
mkdir -p logs
N=$(grep -cv '^\s*\(#\|$\)' "$LIST")
MAXRUN="${HCF_MAXRUN:-20}"
AID=$(sbatch --parsable --array=1-"$N"%"$MAXRUN" "$@" simulation/cluster/run_array.sbatch "$LIST")
BID=$(sbatch --parsable --dependency=afterany:"$AID" simulation/cluster/analyse.sbatch)
echo "$(date -Is) list=$LIST runs=$N array=$AID analyse=$BID" | tee -a logs/submitted.txt
