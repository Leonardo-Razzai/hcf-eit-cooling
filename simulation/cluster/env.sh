# Environment for jobs on the cluster; sourced by the *.sbatch scripts.
# Edit the two variables below once, after installing (see README.md, "Set-up").

# Path of the virtual environment created at installation:
HCF_VENV="${HCF_VENV:-$HOME/venvs/hcf-eit}"
# Modules to load before activating it (leave empty if the system python3 is >= 3.10):
HCF_MODULES="${HCF_MODULES:-}"

if [ -n "$HCF_MODULES" ] && command -v module >/dev/null 2>&1; then
    module load $HCF_MODULES
fi
# shellcheck disable=SC1091
source "$HCF_VENV/bin/activate"

# The 2D Zeeman code is BLAS-bound (dense matrix products): one BLAS thread per allocated core.
NT="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
export MPLBACKEND=Agg
mkdir -p logs
