# Fair PQ Operating Envelopes under Voltage-Unbalance Constraints

Julia code for computing **fair PQ operating envelopes (OEs)** for single-phase PV units in
three-phase unbalanced low-voltage distribution networks, with optional **voltage
unbalance factor (VUF)** limits. Two fairness methods are compared against a
utilitarian benchmark, and the efficiency loss of each is reported as a **Price of
Fairness (PoF)**.

- **M1: sensitivity-weighted fairness.** Every unit moves along a common direction, with step lengths proportional to the inverse of its voltage sensitivity.
- **M2: lexicographic max-min fairness.** Units share a common step. The most constrained units are frozen stage by stage, and the remaining units keep expanding.
- **Utilitarian benchmark.** Each unit gets its own step, and the sum of steps is maximised. This gives the upper bound used for the PoF.

All problems are solved on an exact three-phase AC-OPF (rectangular voltages,
`S = diag(V)·conj(Y·V)`) with [JuMP](https://jump.dev) and
[Ipopt](https://github.com/coin-or/Ipopt).

---

## Repository structure

```
.
├── src/
│   ├── OEFramework.jl        # module entry point (includes everything below)
│   ├── types.jl              # SystemData, VUFConfig, OPFResult, PVUnit, OEBoundary
│   ├── network.jl            # admittance / impedance matrices, XLSX export
│   ├── opf.jl                # AC-OPF model, power balance, VUF constraints
│   ├── envelopes.jl          # M1 (weighted) and M2 (lex max-min) envelopes
│   ├── fairness.jl           # utilitarian benchmark and Price of Fairness
│   ├── reporting.jl          # console result tables
│   ├── plotting.jl           # envelope plots and overlays
│   ├── dashboard.jl          # baseline OPF dashboard (P, Q, voltage, line flows)
│   ├── pipeline.jl           # run_case / run_scenario / run_price_of_fairness
│   └── networks/
│       ├── nw1_13bus.jl      # NW1: modified IEEE 13-node feeder, 400 V
│       └── nw2_116bus.jl     # NW2: IEEE European LV test feeder (116-node), 416 V
├── OE13_6.jl                 # run script, NW1 (13-bus)
├── Final_fairness_55.jl      # run script, NW2 (116-bus)
├── Project.toml              # Julia dependencies
└── README.md
```

---

## Requirements

- Julia **1.10 or later** (developed and tested with 1.11.4).
- Packages: JuMP, Ipopt, MathOptInterface, DataFrames, Plots (GR backend), Colors and XLSX.
  - They are listed in `Project.toml`.
  - Ipopt ships with the MUMPS linear solver through `Ipopt_jll`, so no commercial solver or license is needed.

## Installation

```bash
git clone https://github.com/alireza33zz/OE_VU-F2-v1.0.git
cd OE_VU-F2-v1.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

Run either script from the repository root, with the project environment active:

```bash
julia --project=. OE13_6.jl            # NW1, 13-bus  (a few minutes)
julia --project=. Final_fairness_55.jl # NW2, 116-bus (considerably longer)
```

From the Julia REPL or VS Code you can instead activate the environment (`] activate .`) and run the script.

Each script runs these steps:

1. Builds the network and runs short warm-up calls to trigger JIT compilation.
2. Solves four scenarios: M1 and M2, each without and with VUF limits. For each scenario it:
   - solves the baseline minimum-cost AC-OPF;
   - computes the per-phase PQ envelopes;
   - saves the dashboard and envelope figures.
3. Prints a runtime summary.
4. Runs the utilitarian benchmark for each VUF case and prints the **Price of Fairness** of M1 and M2.
5. Saves overlay figures of all four scenarios to `Graphics/`.

### Main settings (top of each script)

| Setting | NW1 | NW2 | Meaning |
|---|---|---|---|
| `vuf_threshold` | 0.010 | 0.015 | VUF limit ε (fraction) |
| `vuf_buses` | `[11, 12, 13]` | `[40, 60, 80, 96]` | buses where the limit is enforced |
| `Angle_step` | 5° | 5° | direction step (72 directions) |
| `m2_kw` | 3000 / 5000 | 3000 / 5000 | Ipopt iteration caps for the M2 intermediate / final stages |
| `m1_kw` | – | `voltage_ref = :unit` | voltage used in the M1 sensitivity weights |
| `TEST_NUM` | 7 | 4 | tag appended to output folder and file names |

### Output

| Path | Content |
|---|---|
| `opf_fairness_figures/<case>/` | M1: baseline dashboard (`opf_*.pdf`) and envelopes (`pq_OE.pdf`) |
| `opf_fairness_LEX_figures/<case>/` | M2: same figures |
| `Graphics/overlay_M1_M2_states_TEST*.{pdf,svg}` | overlay of all four scenarios |
| `system_matrix.xlsx` | NW1 admittance matrix |
| console | baseline OPF results, runtimes, PoF table |

### Using the module directly

```julia
include("src/OEFramework.jl")
using .OEFramework

sys = build_system_data_13bus()
r1  = run_case(sys; method = :weighted,  n_oe_angles = 36)
r2  = run_case(sys; method = :lexmaxmin, n_oe_angles = 36)
pof = run_price_of_fairness(sys, r1, r2; n_oe_angles = 36)
pof.M1, pof.M2
```

To use your own network, write a function that returns a `SystemData`, following `src/networks/nw1_13bus.jl`.

---

## Test networks

- **NW1:** modified IEEE 13-node feeder, 400 V line-to-line, with 9 single-phase PV units (3 per phase).
- **NW2:** IEEE European LV test feeder reduced to 116 nodes, 416 V line-to-line.

The network data are defined in Julia in `src/networks/`.

## Method summary

The network is a three-phase AC-OPF with:

- a balanced slack at bus 1;
- voltage limits;
- PV apparent-power limits $P_u^2+Q_u^2 \le S_u^{\max 2}$;
- optionally, VUF limits $|V^-|^2 \le \varepsilon^2 |V^+|^2$ at the selected buses.

The baseline $(\hat P_u,\hat Q_u)$ is the minimum-cost OPF solution.

Envelopes are traced **per phase**: units on the other phases stay fixed at their baseline. For each direction $\theta_m = 2\pi m/M$, every unit $u$ on that phase moves along

$$P_u = \hat P_u + \delta_u \cos\theta_m, \qquad Q_u = \hat Q_u + \delta_u \sin\theta_m, \qquad \delta_u \ge 0 .$$

| Method | Step rule | Objective |
|---|---|---|
| M1 | $\delta_u = w_u\,\alpha$, with $w_u \propto 1/\lVert\partial V/\partial S_u\rVert$ normalised per phase | $\max \alpha$ |
| M2 | common $\alpha$ for free units, staged freezing of the unit with the smallest capacity slack | $\max \alpha$ at each stage |
| Utilitarian | individual $\delta_u$ | $\max \sum_u \delta_u$ |

**Price of Fairness** (reported value), where $A_u$ is the area of unit $u$'s envelope polygon:

$$\mathrm{PoF}^k = 1 - \frac{\sum_u A_u^{k}}{\sum_u A_u^{\text{util}}}.$$

The utilitarian problem is a nonconvex NLP, so it is solved to local optimality. It is warm-started from the better of M1 and M2 at each direction. The code warns if the fair solution ever beats the benchmark.

### Known limitations

- **Local optimality.** All problems are nonconvex, so Ipopt returns local optima.
- **Failed directions.** A direction that fails to solve is recorded at the baseline point (zero step) so that the angles stay aligned. Failures are printed to the console.
- **Runtime.** NW2 with 72 directions takes a long time, which is dominated by M2. Reduce `Angle_step` resolution for quick tests.
- **Plot font warning.** The GR backend may print `GKS: could not find font bold.ttf`. The warning is harmless.

---

## Citation

The accompanying paper is under review. A citation will be added here once it is published.

## License

A license has not been added yet. Until one is added, all rights are reserved by the author.

## Contact

Alireza Zabihi, Universidad Politécnica de Madrid.
Please use [GitHub Issues](../../issues) for bug reports and questions.
