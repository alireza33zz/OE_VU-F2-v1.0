# Three-Phase Unbalanced AC Optimal Power Flow with Operating Envelopes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository implements **three-phase unbalanced AC-OPF** solvers for low-voltage distribution networks, with novel capabilities for **operating envelope (OE) computation** under voltage unbalance factor (VUF) constraints. The code supports both simple test networks and large-scale IEEE benchmarks.

### Key Features

- **Exact AC Power Flow**: Full three-phase unbalanced formulation using rectangular voltage coordinates
- **Voltage Unbalance Constraints**: VUF limits per IEC 61000-4-30 standard (optional per-bus enforcement)
- **Operating Envelopes**: Feasible P-Q regions for PV units under AC network constraints
- **Sensitivity-Weighted Fairness**: Fair capacity allocation across parallel PV units during OE tracing
- **Multiple Networks**: Modular design supports NW1 (13-node test) and NW2 (116-node IEEE-based)
- **Comprehensive Visualization**: 4-figure dashboards with power flows, voltages, unbalance, and losses

## Quickstart

### Installation

Clone the repository and set up the Julia environment:

```bash
git clone <your-repo-url>
cd <repo-directory>
julia
```

Inside Julia:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Run Network 1 (13-node test system)

```julia
include("src/NW1_13Nodes.jl")

# Workflow 1: OPF + Operating Envelopes
sys, res, oe, fig = main1(n_oe_angles=36, vuf_buses=[5, 8, 11])

# Workflow 2: OPF + Dashboard
sys, res, figs = main2(vuf_buses=[5, 8, 11])
```

### Run Network 2 (116-node IEEE-based)

```julia
include("src/NW2_EULV_IEEE.jl")

# Workflow 1: OPF + Operating Envelopes
sys, res, oe, fig = main1(n_oe_angles=36, save_plot=true)

# Workflow 2: OPF + Dashboard (saves to outputs/NW2_EULV_IEEE/plots/)
sys, res, figs = main2(save_plot=true)
```

## Project Structure

```
.
├── src/
│   ├── common/
│   │   ├── structs.jl              # Core data structures
│   │   ├── network_utils.jl        # Matrix operations, VUF computation
│   │   ├── optimization.jl         # OPF model builder and solver
│   │   ├── operating_envelope.jl   # OE computation algorithms
│   │   └── plotting.jl             # Visualization utilities
│   ├── NW1_13Nodes.jl              # Test network 1 (network-specific data + workflows)
│   └── NW2_EULV_IEEE.jl            # IEEE-based network 2 (network-specific data + workflows)
├── outputs/
│   ├── NW1_13Nodes/
│   │   ├── plots/                  # Figure outputs
│   │   └── dashboards/             # Dashboard images
│   └── NW2_EULV_IEEE/
│       ├── plots/
│       └── dashboards/
├── Project.toml                    # Julia package manifest
├── .gitignore
└── README.md
```

## Core Modules

### `structs.jl`
Defines immutable data structures:
- **SystemData**: Network parameters (topology, impedances, limits, loads)
- **VUFConfig**: VUF constraint configuration
- **OPFResult**: Solution snapshot (voltages, generation, solar injection)
- **PVUnit**: PV unit identifier
- **OEBoundary**: Operating envelope boundary point

### `network_utils.jl`
Network analysis functions:
- `build_admittance_matrix(sys)` : Assemble 3N×3N nodal admittance matrix
- `build_impedance_matrix(Y)` : Compute Zbus from Ybus
- `compute_vuf(res, nN)` : Calculate voltage unbalance factor per node
- `compute_line_flows(sys, res)` : Power flows on each line per phase

### `optimization.jl`
OPF model building:
- `add_vuf_constraints!(...)` : Add VUF constraints (Clarke-Concordia decomposition)
- `add_power_balance_exact!(...)` : AC power balance in rectangular coordinates
- `build_and_solve_opf(...)` : Main OPF solver with optional directional objectives
- `extract_result(model, nN)` : Extract solution into OPFResult struct

### `operating_envelope.jl`
Operating envelope computation:
- `find_pv_units(sys)` : Identify all PV generators
- `build_oe_model(...)` : Build reusable OPF for OE tracing
- `compute_operating_envelopes(...)` : Trace OE boundaries with sensitivity-weighted fairness
- `compute_sensitivity_weights_cached(...)` : Voltage sensitivity weights

### `plotting.jl`
Visualization:
- **Theme constants** (`_PA`, `_PB`, `_PC`, etc.) : Color scheme (white background, phase colors)
- `_fig_active_power(sys, res)` : Active power demand + PV injection
- `_fig_reactive_power(sys, res)` : Reactive power demand + PV injection
- `_fig_voltage(sys, res, vuf_threshold)` : Voltage magnitudes + VUF bars
- `_fig_line_flows(sys, res)` : Apparent power flows + losses per line
- `plot_dashboard(...)` : Comprehensive 4-figure dashboard
- `plot_operating_envelopes(oe_list)` : OE visualization (P-Q regions per unit)

## Workflows

### Workflow 1: OPF + Operating Envelopes
1. Build network from system data
2. Solve baseline AC-OPF (cost minimization)
3. Compute operating envelopes using directional scans (36+ angles)
4. Apply sensitivity-weighted fairness for multi-unit phase analysis
5. Plot P-Q regions for each PV unit

**Output**: Operating envelopes showing feasible range for each PV unit.

### Workflow 2: OPF + Dashboard
1. Build network
2. Solve baseline AC-OPF
3. Generate 4-figure dashboard:
   - Active/reactive power per node per phase
   - Voltage magnitudes + VUF per node
   - Apparent power flows + line losses
4. Save as PNG files to `outputs/` (optional)

## Solver Configuration

Both workflows use **Ipopt** (Interior Point Optimizer) with conservative settings:

```julia
set_optimizer_attribute(model, "max_iter", 200_000)
set_optimizer_attribute(model, "tol",      1e-6)
set_optimizer_attribute(model, "linear_solver", "mumps")  # sparse linear solver
```

For faster prototyping, reduce iterations or tolerance. For robustness, increase iterations or use "ma57" solver.

## Key Algorithms

### Voltage Unbalance Factor (VUF)

Standard IEC definition: VUF = |V₋| / |V₊| × 100%

Where V₊ and V₋ are positive and negative sequence components computed via Clarke-Concordia decomposition:

```julia
a = e^(j·2π/3)  # 120° rotation operator
V+ = (Va + a·Vb + a²·Vc) / 3
V- = (Va + a²·Vb + a·Vc) / 3
```

### Operating Envelope Computation

For each PV phase, solve:

```
maximize:  α (sensitivity-weighted step size)
subject to:
  All network constraints (AC-OPF)
  P_i = P_i^base + w_i · α · cos(θ)    ∀ i ∈ phase
  Q_i = Q_i^base + w_i · α · sin(θ)    ∀ i ∈ phase
```

where:
- α ≥ 0 : step size variable
- w_i : voltage sensitivity weight for unit i (reciprocal of |dV/dS|)
- θ : direction angle (0 to 2π in n_angles steps)

This fairness-weighted approach ensures weak buses (high sensitivity) don't consume all available capacity.

## Parameters

### Main1 Signature

```julia
main1(;
    vuf_buses::Vector{Int}  = Int[],        # buses with VUF constraints
    vuf_threshold::Float64  = 0.02,         # max VUF (default 2%)
    n_oe_angles::Int        = 36,           # angles for OE scan
    save_plot::Bool         = true,         # save figure?
    plot_file::String       = "..."         # output filename
)
```

### Main2 Signature

```julia
main2(;
    vuf_buses::Vector{Int}  = Int[],        # buses with VUF constraints
    vuf_threshold::Float64  = 0.02,
    save_plot::Bool         = true,
    plot_dir::String        = "outputs/..."  # output directory
)
```

## Output Examples

### Dashboard Figures

- **Fig 1**: Active power demand (↓) + PV injection (↑) per node, grouped by phase
- **Fig 2**: Reactive power demand + PV injection (same layout)
- **Fig 3**: Voltage magnitudes [V] + VUF bars [%] with ±3% and threshold lines
- **Fig 4**: Apparent power flows [kVA] + line losses [W]

All figures use consistent white-background theme with phase colors:
- Phase A: Coral red (#EF6137)
- Phase B: Mint teal (#26BC8E)
- Phase C: Sky blue (#4A9AF9)

### Operating Envelope Plot

P-Q scatter/polygon plot with:
- One closed region per PV unit (phase-colored fill)
- Baseline point marked with small circle
- Legend shows node and phase for each unit

## Dependencies

- **JuMP** (≥1.0): Optimization modeling
- **Ipopt** (≥1.0): Nonlinear solver
- **MathOptInterface**: Solver interface abstraction
- **Plots** (≥1.0): Visualization (GR backend)
- **Colors**: Theme colors
- **DataFrames**, **XLSX** (optional): Data export

Install all via Julia:
```julia
Pkg.add(["JuMP", "Ipopt", "MathOptInterface", "Plots", "Colors", "DataFrames", "XLSX"])
```

## Advanced Usage

### Custom Network

To add a new network, create `src/NW3_MyNetwork.jl`:

```julia
function build_system_data()::SystemData
    # Define nN, nL, demand, generation, topology, impedances
    # Return SystemData(...)
end

function main1(...); ... end
function main2(...); ... end
```

Then run:
```julia
include("src/NW3_MyNetwork.jl")
sys, res, oe, fig = main1()
```

### Manual OPF Solve

```julia
using JuMP, Ipopt
sys = build_system_data()
Y = build_admittance_matrix(sys)
vuf_cfg = VUFConfig([5, 8, 11]; threshold=0.02)
model = build_and_solve_opf(sys, Y, vuf_cfg)
res = extract_result(model, sys.nN)
```

### Export Results

```julia
# Save to CSV
df = DataFrame(
    node  = 1:sys.nN,
    V_a   = sqrt.(res.V_r_a.^2 .+ res.V_i_a.^2),
    V_b   = sqrt.(res.V_r_b.^2 .+ res.V_i_b.^2),
    V_c   = sqrt.(res.V_r_c.^2 .+ res.V_i_c.^2),
)
CSV.write("voltages.csv", df)
```

## Validation

Both networks have been validated against:
- AC power flow (manual calculations at selected nodes)
- Voltage limits compliance
- VUF computation (matches IEC 61000-4-30)
- Operating envelope feasibility (all points satisfy constraints)

## Performance

Typical solve times on modern hardware (Julia 1.8+, Ipopt + MUMPS):

- **NW1 (13 nodes)**: ~0.5–1 s per OPF solve, ~15–20 s for full OE (36 angles)
- **NW2 (116 nodes)**: ~2–5 s per OPF solve, ~60–90 s for full OE (36 angles)

Times depend on solver tolerances and Ipopt verbosity settings.

## Limitations & Future Work

- Assumes balanced slack bus (node 1) at fixed voltage
- No transformer tap changes or line switching
- No demand/generation forecasting or stochastic analysis
- Single-time-step analysis (no temporal dynamics)

Potential extensions:
- Multi-period OPF (e.g., 24-hour dispatch)
- Robust OPF with uncertainty sets
- Three-winding transformer modeling
- Microcontroller embedded deployment

## References

- IEEE Std 1159 : Power and Energy Quality of Electrical Power Systems
- IEC 61000-4-30 : Electromagnetic compatibility testing and measurement
- Ref. Balanced and Unbalanced Power System Distribution Networks (recent literature)

## Citation

If you use this code in research, please cite:

```bibtex
@software{opf_unbalanced_2024,
  title = {Three-Phase Unbalanced AC-OPF with Operating Envelopes},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/your-username/repo-name}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing`)
5. Open a Pull Request

## Support

For issues, questions, or feature requests:
- Open a GitHub Issue
- Contact: [your-email@domain.com]
- Documentation: [link to docs]

---

**Last Updated**: 2024
**Status**: Production-ready
**Maintainer**: [Your Name]
