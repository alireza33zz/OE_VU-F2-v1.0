# Operating Envelopes — 3-Phase Unbalanced AC-OPF

Lexicographic max-min and weighted-fairness **PQ Operating Envelopes** for
3-phase unbalanced distribution networks with Voltage Unbalance Factor (VUF)
constraints.

## Repository structure

```
oe_repo/
├── src/
│   ├── types.jl                # Structs: SystemData, VUFConfig, OPFResult,
│   │                           #          PVUnit, OEBoundary
│   ├── admittance.jl           # build_admittance_matrix, build_impedance_matrix
│   ├── opf_model.jl            # add_vuf_constraints!, add_power_balance_exact!,
│   │                           # build_and_solve_opf (V_base kwarg),
│   │                           # build_oe_model (V_base kwarg),
│   │                           # set_directional_objective!, extract_result
│   ├── postprocess.jl          # compute_vuf, compute_line_flows
│   ├── operating_envelopes.jl  # find_pv_units, angular_diff, local_optimal_angle,
│   │                           # compute_sensitivity_weights_cached,
│   │                           # compute_operating_envelopes (fairness),
│   │                           # compute_operating_envelopes_lexmaxmin
│   └── plotting.jl             # print_results, plot_dashboard, _fig_*,
│                               # plot_operating_envelopes_overlay_legend,
│                               # plot_operating_envelopes_overlay
├── networks/
│   ├── nw1_system_data.jl      # 13-node, 400 V (IEEE modified feeder)
│   └── nw2_system_data.jl      # EULV 116-node, 416 V (IEEE European LV)
├── outputs/                    # Auto-created; all PDFs written here (git-ignored)
├── main.jl                     # Single entry point — runs both test cases
├── Project.toml                # Julia package dependencies
└── .gitignore
```

## How to run

### From the Julia REPL (recommended)

```julia
# Activate the project environment first
using Pkg; Pkg.activate(".")
Pkg.instantiate()          # first time only — downloads all dependencies

include("main.jl")
```

### From the terminal

```bash
julia --project=. main.jl
```

Figures are exported as PDF to `outputs/NW1/` and `outputs/NW2/`.

## Configuring a study

Edit the **Study parameters** block near the bottom of `main.jl`:

```julia
VUF_THRESHOLD = 0.010   # IEC 61000-3-13 limit (0.01 = 1%, 0.02 = 2%)
ANGLE_STEP    = 5       # Tracing resolution in degrees (5 -> 72 directions)
```

Adjust each `run_case(...)` call to change which buses have VUF enforced:

```julia
run_case("NW1", build_system_data_nw1(), NW1_V_BASE;
    vuf_buses     = [11, 12, 13],   # buses where VUF <= threshold
    vuf_threshold = VUF_THRESHOLD,
    n_oe_angles   = N_ANGLES)
```

## Adding a new network

1. Create `networks/nwX_system_data.jl` with a `build_system_data_nwX()` function
   returning a `SystemData` and a `NWX_V_BASE` constant.
2. Add one `run_case(...)` call in `main.jl` — no changes to `src/` needed.

## Key design notes

| Topic | Detail |
|---|---|
| `V_base` kwarg | `build_and_solve_opf`, `build_oe_model`, and both OE functions accept `V_base::Float64` (default `400/sqrt(3)`). Pass `NW1_V_BASE` or `NW2_V_BASE` from the network file. |
| Sensitivity indexing | `compute_operating_envelopes_lexmaxmin` uses the corrected per-unit sensitivity: `z_col ./ conj(V_base_unit)` (NW2 version). |
| Plot layout | `plot_operating_envelopes_overlay_legend` uses the compact IEEE two-column layout (NW2 version). |
| Solver tolerances | `max_iter` and `tol` are keyword arguments — NW1 uses tighter tolerances by default (`max_iter=200_000, tol=1e-6`); NW2 uses looser ones (`max_iter=10_000, tol=1e-4`). Override in `run_case` if needed. |

## Dependencies

All managed via `Project.toml`. Key packages:

- [JuMP.jl](https://jump.dev/) — algebraic modelling
- [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) — nonlinear interior-point solver
- [Plots.jl](https://docs.juliaplots.org/) — all figures (GR backend)
- [XLSX.jl](https://felipenoris.github.io/XLSX.jl/) — admittance matrix export (NW1)
