#!/bin/bash
# From your computer: status of the cluster jobs and fetch of the 2D Zeeman results.
#     bash simulation/cluster/remote.sh status     # queue, number of finished runs, analysis done?
#     bash simulation/cluster/remote.sh fetch      # rsync results/zeeman2d (+ zeeman2d.dat) back, rerun the plots locally
#     bash simulation/cluster/remote.sh summary    # print the cluster-side summary.txt
# Needs a non-interactive ssh to the frontend: an SSH key accepted by the bastion, or an open
# ControlMaster connection (README.md, "Automatic checks"). Settings:
HOST="${HCF_HOST:-ophfe1}"
REMOTE="${HCF_REMOTE:-hcf-eit-cooling}"          # repository path on the cluster, relative to the home
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SSH="ssh -o BatchMode=yes -o ConnectTimeout=15"
case "${1:-status}" in
  status)
    $SSH "$HOST" "cd $REMOTE && echo '--- queue' && squeue -u \$USER -o '%.12i %.20j %.8T %.10M %.6D %R' && \
      echo '--- finished runs:' \$(ls simulation/results/zeeman2d/runs/*.json 2>/dev/null | wc -l) && \
      echo '--- last submissions' && tail -3 logs/submitted.txt 2>/dev/null; \
      if [ -f simulation/results/zeeman2d/summary.txt ]; then echo '--- summary.txt from' \$(date -r simulation/results/zeeman2d/summary.txt); fi"
    ;;
  fetch)
    rsync -az -e "$SSH" "$HOST:$REMOTE/simulation/results/zeeman2d/" "$ROOT/simulation/results/zeeman2d/"
    rsync -az -e "$SSH" "$HOST:$REMOTE/simulation/results/zeeman2d.dat" "$ROOT/simulation/results/" 2>/dev/null || true
    cd "$ROOT/simulation/plotting" && "$ROOT/.venv/bin/python" zeeman2d_plots.py
    ;;
  summary)
    $SSH "$HOST" "cat $REMOTE/simulation/results/zeeman2d/summary.txt"
    ;;
  *) echo "usage: remote.sh status|fetch|summary"; exit 1 ;;
esac
