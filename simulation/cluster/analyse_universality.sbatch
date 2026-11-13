#!/bin/bash
# Collect and plot the universality runs (studies/universality.py, plotting/universality_plots.py).
#     sbatch --dependency=afterany:<jobs> simulation/cluster/analyse_universality.sbatch
#SBATCH --job-name=hcf_univ_analyse
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source simulation/cluster/env.sh
cd simulation/studies
python universality.py collect
python ../plotting/universality_plots.py
