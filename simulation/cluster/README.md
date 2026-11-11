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

## Set-up (once)

```bash
# on the cluster, in your home or work directory
git clone <repository URL> hcf-eit-cooling        # or: rsync -a --exclude .venv ./ cluster:hcf-eit-cooling/
cd hcf-eit-cooling
python3 -m venv ~/venvs/hcf-eit                    # Python >= 3.10 (load a module first if needed)
source ~/venvs/hcf-eit/bin/activate
pip install -r simulation/requirements.txt
```

Then edit `simulation/cluster/env.sh` if the venv path differs, or if a module must be loaded (`HCF_MODULES`). Submit every job **from the repository root**: logs go to `logs/`, results to `simulation/results/zeeman2d/runs/`.

## 1. Check the installation (a few minutes, debug queue)

```bash
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py check
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py check_field
sbatch --qos=debug --time=00:15:00 simulation/cluster/run_single.sbatch zeeman2d.py smoke
```

- `check` and `check_field` compare the 2D code (with N_y = 1) with the validated one-axis Zeeman solver: at zero field, and at μ_B B/ħ = 40 ν_r. Both must print `PASS`.
- `smoke` builds all four schemes and takes a few steps.

Then a two-line test scan at N = 4 (≈ 5 min per line):

```bash
sbatch --array=1-2 --time=00:30:00 simulation/cluster/run_array.sbatch simulation/cluster/scans/zeeman2d_test.txt
```

## 2. Main scan (216 runs at N = 6)

```bash
LIST=simulation/cluster/scans/zeeman2d_main.txt
N=$(grep -cv '^\s*\(#\|$\)' $LIST)
sbatch --array=1-$N%20 simulation/cluster/run_array.sbatch $LIST
```

The scan covers, for cases A and B at s = 0.5 and 1:
- **TE₀₁ (both field directions):** two-photon detunings d₁, d₂ ∈ {0, 0.4, 0.8}, with and without repumper.
- **`lp11_hf`:** d ∈ {0.3, 0.6, 0.9}, repumper 0.03 or 0.1, stage length τ = 16 or 32.
- **`lp11_zl`:** d ∈ {−0.4, 0, 0.4}, repumper 0.3 or 1.

Runs whose JSON already exists are skipped, so resubmit the same command after timeouts. To regenerate or change the list: `cd simulation/studies && python zeeman2d.py list --grid main [--N 6]`, or edit `main_grid()`.

## 3. Convergence in N (long queue)

```bash
LIST=simulation/cluster/scans/zeeman2d_convergence.txt   # te01_bz and lp11_hf, A and B, N = 6, 8, 10
N=$(grep -cv '^\s*\(#\|$\)' $LIST)
sbatch --array=1-$N --qos=long --time=72:00:00 --mem=8G --cpus-per-task=8 simulation/cluster/run_array.sbatch $LIST
```

## 4. Collect and bring back

```bash
cd simulation/studies && python zeeman2d.py collect      # -> simulation/results/zeeman2d.dat
# locally:
rsync -av cluster:hcf-eit-cooling/simulation/results/zeeman2d/ simulation/results/zeeman2d/
```

**Each run writes a JSON file containing:**
- the parameters;
- P₀₀ (mean over the last switching period, or over 50/ν_r);
- W₀₀ = P₀₀ss / ∫(P₀₀ss − P₀₀)dt from |1,1⟩, in ν_r and in s⁻¹;
- `drift`: the change of P₀₀ over the last window, which must be ≪ 10⁻³, otherwise increase `--tmax`;
- n_x, n_y, the sublevel populations;
- `frame_mismatch`: must be ~10⁻¹²;
- P₀₀(t), down-sampled.

The internal state starts uniform over F=2 (`--init target` starts in the target sublevel).

## Single runs

```bash
sbatch simulation/cluster/run_single.sbatch zeeman2d.py run --scheme te01_bz --case B --s 1 --d1 0.4 --d2 0.4 --N 6
```

`python zeeman2d.py run -h` lists all options (Δ/Γ, field `--zeeman` in units of ν_r, coupling scale `--wc`, `--tmax`, `--tau`, ...).

## Notes

- **Field strength.** The field enters as Z = μ_B B/ħ in units of ν_r. The default Z = 40 is 0.29 MHz in case A (B ≈ 0.2 G) and 1.2 MHz in case B (B ≈ 0.9 G). Different Raman resonances are Z/2 apart, and the code merges paths closer than `--tol` (default 2 ν_r). Keep Z/2 ≫ tol ≫ |d|.
- **Results directory.** `simulation/results/zeeman2d/runs/` is created by the runs. Commit the collected `.dat` file (and the JSONs, if wanted) when the results are final.
- **Threads.** Set by `--cpus-per-task` through `env.sh`. Do not run several jobs per core.
