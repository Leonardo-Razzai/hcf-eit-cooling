# Running the 2D Zeeman simulations on a Slurm cluster

`simulation/studies/zeeman2d.py` (module `simulation/eitsim/zeeman2d.py`) solves the eliminated master equation on
(8 ground sublevels of ⁸⁷Rb) × (n_x, n_y). It compares four schemes:

| Scheme | What it is |
|---|---|
| `lp11_hf` | LP₁₁ orientation switched; hyperfine Λ, D1 F′=1, π probe on F=2, ⊥ coupling, repumpers (best one-axis scheme) |
| `lp11_zl` | LP₁₁ switched; Zeeman Λ inside F=2 (σ⁻ probe, σ⁺ coupling, B ∥ fibre, F=1 repumper) |
| `te01_bz` | TE₀₁ probe, B ∥ fibre, σ⁺ coupling on F=1, F′=2, target \|2,1⟩, two probe frequencies |
| `te01_bx` | TE₀₁ probe, B transverse, π coupling on F=1, F′=2, target \|2,1⟩, two probe frequencies |

Physics and level schemes are in the background report, Sec. "One vector mode for both radial axes".

**Cost.** Dense RK4 on a density matrix of dimension d = 8 N², limited by BLAS speed. Measured on 4 threads (laptop), tmax = 1500/ν_r:

| N per axis | d | Time | Memory |
|---|---|---|---|
| 4 | 128 | ≈ 7 min | < 0.5 GB |
| 6 | 288 | ≈ 1 h | < 1 GB |
| 8 | 512 | ≈ 6 h | ≈ 1–2 GB |
| 10 | 800 | ≈ 1 day (estimate) | ≈ 2–4 GB |

`python zeeman2d.py cost` prints flops and operator sizes. One run uses one node share (`--cpus-per-task`); more cores help up to ~8.

## DIFA OPH cluster: rules and access

These notes are taken from `rb_lattice_dsmc/README.md`, which follows the DIFA wiki (startup guide, access, storage and jobs pages, October 2026).

- **Frontends** (`ophfe1`, `ophfe2`) are for logging in, editing, installing and submitting only. Never run simulations or the test suite there (1 GB RAM per user); everything goes through `sbatch`.
- **Login** through the bastion. In `~/.ssh/config` (in WSL: `~/.ssh` of Ubuntu):
  ```
  Host bastion-nav
      HostName bastion-nav.difa.unibo.it
      User name.surname00
  Host ophfe1 ophfe2
      User name.surname00
      ProxyJump bastion-nav
  ```
  then `ssh ophfe1`.
- **Queues:** `--qos=debug` (≤ 15 min, not billed) for the checks, `--qos=long` (≤ 72 h, low priority) for N ≥ 8. Ask your sector contact which node type to use (`--constraint=blade` or `matrix`) and whether `--account` is needed. Then uncomment those `##SBATCH` lines in `run_single.sbatch` and `run_array.sbatch`.
- **Support:** access and accounting through your sector contact; technical problems through `difa.csi@unibo.it`.

## Set-up (once, on the frontend)

```bash
git clone https://github.com/Leonardo-Razzai/hcf-eit-cooling.git ~/hcf-eit-cooling   # private repo: token or SSH key
cd ~/hcf-eit-cooling
python3 --version                    # need >= 3.10 (numpy 2)
module avail 2>&1 | grep -i python   # if the system python is older, look for a module
# module load <python module>        # if you use one, also put it in HCF_MODULES in simulation/cluster/env.sh
python3 -m venv ~/venvs/hcf-eit
source ~/venvs/hcf-eit/bin/activate
python -m pip install --upgrade pip
python -m pip install -r simulation/requirements.txt
mkdir -p logs                        # Slurm writes the job logs here: it must exist BEFORE the first sbatch
```

- **No package installation is needed:** the scripts put `simulation/` on `sys.path` themselves (`studies/_common.py`), and the job scripts run them from `simulation/studies`.
- **Keep the venv in `/home`.** It is backed up (50 GB quota); `/scratch` is not.
- **Without a Python ≥ 3.10 module,** install Miniforge in your home (`bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3`) and create the venv with `~/miniforge3/bin/python3 -m venv ~/venvs/hcf-eit`.
- **If the venv is elsewhere,** set `HCF_VENV` at the top of `simulation/cluster/env.sh`.
- **Updates later:** `git pull` (or copy `simulation/` with `scp -r`).

Submit every job **from the repository root** (`cd ~/hcf-eit-cooling`): logs go to `logs/`, results to `simulation/results/zeeman2d/runs/`.

## 1. Check the installation (debug queue, a few minutes)

```bash
# regression tests on a compute node (~30 s; must end with "35 passed")
sbatch --qos=debug --time=00:15:00 --output=logs/tests_%j.out \
       --wrap "source simulation/cluster/env.sh && python -m pytest simulation/tests -q"
# 2D code against the one-axis Zeeman solver, zero and finite field: both must print PASS
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py check
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py check_field
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py smoke
squeue -u $USER                      # or: slurmtop
cat logs/hcf2d_<jobid>.out
```

Then a two-line test scan at N = 4, case A (≈ 5 min per line):

```bash
sbatch --array=1-2 --time=00:30:00 simulation/cluster/run_array.sbatch simulation/cluster/scans/zeeman2d_test.txt
```

## 2. Main scan (108 runs at N = 6, case A)

```bash
bash simulation/cluster/submit_scan.sh simulation/cluster/scans/zeeman2d_main.txt
```

This submits the array, at most 20 tasks at a time (`HCF_MAXRUN` changes this). It then submits `analyse.sbatch` with `--dependency=afterany:<array>`. As soon as the last run ends, that job:
- collects the runs into `simulation/results/zeeman2d.dat`;
- writes `simulation/results/zeeman2d/summary.txt`, `overview`, `detunings`, `convergence` and `traces` (pdf/png).

The job ids are appended to `logs/submitted.txt`. To get an e-mail when the analysis is ready, uncomment the two `--mail` lines in `analyse.sbatch`.

The scan covers case A (a = 15 µm) only, at s = 0.5 and 1:
- **TE₀₁ (both field directions):** two-photon detunings d₁, d₂ ∈ {0, 0.4, 0.8}, with and without repumper.
- **`lp11_hf`:** d ∈ {0.3, 0.6, 0.9}, repumper 0.03 or 0.1, stage length τ = 16 or 32.
- **`lp11_zl`:** d ∈ {−0.4, 0, 0.4}, repumper 0.3 or 1.

Runs whose JSON already exists are skipped, so resubmit the same command after timeouts. To regenerate or change the list: `cd simulation/studies && python zeeman2d.py list --grid main [--N 6]`, or edit `main_grid()`.

By hand, without the helper: `sbatch --array=1-$N%20 simulation/cluster/run_array.sbatch $LIST`, then `sbatch simulation/cluster/analyse.sbatch`.

## 3. Convergence in N (long queue)

```bash
bash simulation/cluster/submit_scan.sh simulation/cluster/scans/zeeman2d_convergence.txt \
     --qos=long --time=72:00:00 --mem=8G --cpus-per-task=8      # te01_bz and lp11_hf, case A, N = 6, 8, 10
```

## 4. Results

The analysis job writes `summary.txt`. It contains:
- the flagged runs: |drift| ≥ 10⁻³ (increase `--tmax`) or frame mismatch;
- the best run per scheme with P₀₀ ≥ 0.85: P₀₀ and W₀₀ in ν_r and in s⁻¹;
- the closed-Λ two-axis reference of the minimal model;
- convergence in N.

The figures are in the same folder. By hand: `cd simulation/studies && python zeeman2d.py collect [--out DIR] && python ../plotting/zeeman2d_plots.py [--runs DIR]`.

**Each run writes a JSON file containing:**
- the parameters;
- P₀₀ (mean over the last switching period, or over 50/ν_r);
- W₀₀ = P₀₀ss / ∫(P₀₀ss − P₀₀)dt from |1,1⟩, in ν_r and in s⁻¹;
- `drift`, n_x, n_y, the sublevel populations, `frame_mismatch` (must be ~10⁻¹²);
- P₀₀(t), down-sampled.

The internal state starts uniform over F=2 (`--init target` starts in the target sublevel).

## 5. Automatic checks from your computer

`simulation/cluster/remote.sh` checks the queue and fetches results:

| Command | What it does |
|---|---|
| `bash simulation/cluster/remote.sh status` | queue, number of finished runs, whether `summary.txt` exists |
| `bash simulation/cluster/remote.sh summary` | prints the cluster-side summary |
| `bash simulation/cluster/remote.sh fetch` | rsyncs `results/zeeman2d/` back and reruns the plots locally |

It needs ssh without a password prompt. The bastion asks for the UniBo password, so keep one authenticated connection open with OpenSSH connection sharing. Add to `~/.ssh/config` (in WSL):

```
Host bastion-nav ophfe1 ophfe2
    ControlMaster auto
    ControlPath ~/.ssh/cm-%C
    ControlPersist 12h
```

Then log in once with `ssh ophfe1` (password twice) and exit. For the next 12 h, `remote.sh`, and Claude Code in this repository, can reach the cluster without asking again. Claude Code can then poll `remote.sh status` while a session is open, fetch the results when `summary.txt` appears, and analyse them. Settings: `HCF_HOST` (default `ophfe1`), `HCF_REMOTE` (repository path in the cluster home, default `hcf-eit-cooling`).

## Single runs

```bash
sbatch simulation/cluster/run_single.sbatch zeeman2d.py run --scheme te01_bz --case A --s 1 --d1 0.4 --d2 0.4 --N 6
```

`python zeeman2d.py run -h` lists all options (Δ/Γ, field `--zeeman` in units of ν_r, coupling scale `--wc`, `--tmax`, `--tau`, ...).

## Notes

- **Field strength.** The field enters as Z = μ_B B/ħ in units of ν_r. The default Z = 40 is μ_B B/h = 0.29 MHz in case A (B ≈ 0.2 G). Different Raman resonances are Z/2 apart, and the code merges paths closer than `--tol` (default 2 ν_r). Keep Z/2 ≫ tol ≫ |d|.
- **Results on `/scratch`:** append `--out /scratch/<sector>/name.surname00/zeeman2d` to the `run_array.sbatch` line (extra arguments are passed to every run), and also to `collect`. `/scratch` is fast but not backed up, and files older than 40 days are deleted.
- **Results directory.** `simulation/results/zeeman2d/runs/` is created by the runs. Commit the collected `.dat` file (and the JSONs, if wanted) when the results are final.
- **Threads.** Set by `--cpus-per-task` through `env.sh`. Do not run several jobs per core.
