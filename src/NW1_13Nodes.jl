# ============================================================
# NETWORK 1: 13-NODE TEST SYSTEM
# ============================================================
# Three-phase unbalanced AC-OPF with dynamic PQ operating envelopes
# and voltage unbalance factor (VUF) constraints.

using JuMP, Ipopt, DataFrames, Printf, Colors, Plots
import MathOptInterface as MOI
gr()

# Load shared modules
include("../common/structs.jl")
include("../common/network_utils.jl")
include("../common/optimization.jl")
include("../common/operating_envelope.jl")
include("../common/plotting.jl")

# ============================================================
# NETWORK-SPECIFIC: SYSTEM PARAMETERS
# ============================================================

"""
    build_system_data() -> SystemData

Build the 13-node test network with three-phase unbalanced loads.

This network features:
- 13 nodes, 12 distribution lines
- Highly unbalanced demand: Phase A heavily loaded, Phase C lightly loaded
- PV generators at nodes 5, 8, 11 (one per phase)
- Voltage limits: ±10% of nominal (217–254 V LN)

Returns a SystemData struct with all network parameters.
"""
function build_system_data()::SystemData
    nN = 13;  nL = 12

    # Demand scaling factor (kept at 1.0 for target VUF range)
    Demand_scaling = 1.0

    # Three-phase active power demand [W] — highly unbalanced
    # Heavy on Phase A, moderate on Phase B, light on Phase C
    P_Demand_a = Float64.([ 0,    0,    0,    0, 3000,    0,    0, 3000,    0,    0, 3000,    0,    0]) * Demand_scaling
    P_Demand_b = Float64.([ 0,    0,    0,    0,    0, 2000,    0,    0, 2000,    0,    0, 2000,    0]) * Demand_scaling
    P_Demand_c = Float64.([ 0,    0,    0,    0,    0,    0, 1000,    0,    0, 1000,    0,    0, 1000]) * Demand_scaling

    # Reactive power demand [var] — proportional to active (pf ≈ 0.97)
    Q_Demand_a = P_Demand_a * 0.25
    Q_Demand_b = P_Demand_b * 0.25
    Q_Demand_c = P_Demand_c * 0.25

    # Solar PV capacity per phase [VA] — symmetric placement
    S_Sol_Gen_ub_a = Float64.([ 0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500,    0,    0])
    S_Sol_Gen_ub_b = Float64.([ 0,    0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500,    0])
    S_Sol_Gen_ub_c = Float64.([ 0,    0,    0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500])

    # Generation parameters (slack bus at node 1)
    P_Cost_G   = [1.0; zeros(nN-1)]
    P_Gen_lb   = [-3_000_000.0; zeros(nN-1)]
    P_Gen_ub   = [ 3_000_000.0; zeros(nN-1)]
    Q_Gen_lb   = [-500_000.0; zeros(nN-1)]
    Q_Gen_ub   = [ 500_000.0; zeros(nN-1)]
    Gen_Status = [1; zeros(Int,nN-1)]

    # Voltage limits [V]: ±10% of nominal (230 V LN)
    V_base = 400.0 / sqrt(3)
    V_lb   = fill(0.94 * V_base, nN)
    V_ub   = fill(1.10 * V_base, nN)

    # Network topology: from-bus and to-bus for each line
    fbus = [1, 2, 2, 2, 2, 3, 3, 3, 3, 4,  4,  4]
    tbus = [2, 5, 6, 7, 3, 8, 9, 10, 4, 11, 12, 13]

    # Line lengths [km] (scaled for realistic impedances)
    lengths = [0.6, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] .* 0.5

    # Per-unit impedance matrix [Ω/km] (symmetric 3×3 for balanced cable)
    Zij_per_km = ComplexF64[
        0.799+0.075im  0.353+0.004im  0.353+0.004im;
        0.353+0.004im  0.799+0.075im  0.353+0.004im;
        0.353+0.004im  0.353+0.004im  0.799+0.075im]

    # Per-branch impedance matrices: scale by line length
    Zij = [lengths[k] .* Zij_per_km for k in 1:length(fbus)]

    # Shunt admittance (typically zero for distribution feeders)
    Ysh = zeros(ComplexF64, 3, 3)

    return SystemData(nN, nL,
        P_Demand_a, P_Demand_b, P_Demand_c,
        Q_Demand_a, Q_Demand_b, Q_Demand_c,
        P_Cost_G, P_Gen_lb, P_Gen_ub, Q_Gen_lb, Q_Gen_ub, Gen_Status,
        S_Sol_Gen_ub_a, S_Sol_Gen_ub_b, S_Sol_Gen_ub_c,
        V_lb, V_ub, fbus, tbus, Zij, Ysh)
end

# ============================================================
# HELPER: PRINT RESULTS
# ============================================================

"""
    print_results(sys, res, vuf_buses)

Print a formatted summary of OPF solution including voltages, generation,
and VUF at selected buses.
"""
function print_results(sys::SystemData, res::OPFResult, vuf_buses::Vector{Int}=Int[])
    println("\n" * "="^70)
    println("OPF SOLUTION SUMMARY")
    println("="^70)
    
    if res.status != MOI.OPTIMAL && res.status != MOI.LOCALLY_SOLVED
        println("⚠  Solver status: $( res.status)")
        return
    end
    
    println("✓ Status: $(res.status)")
    println("  Objective value: $(round(res.obj; digits=2))")
    
    println("\n── Voltage Magnitudes [V] ─────────────────────────────────────")
    _mag(r, i) = sqrt(r^2 + i^2)
    for i in 1:min(sys.nN, 5)
        Va = _mag(res.V_r_a[i], res.V_i_a[i])
        Vb = _mag(res.V_r_b[i], res.V_i_b[i])
        Vc = _mag(res.V_r_c[i], res.V_i_c[i])
        @printf("  Node %2d:  Va=%.1f  Vb=%.1f  Vc=%.1f\n", i, Va, Vb, Vc)
    end
    
    println("\n── Solar Injection [kW] ──────────────────────────────────────")
    for i in 1:sys.nN
        if sys.S_Sol_Gen_ub_a[i] > 0 || sys.S_Sol_Gen_ub_b[i] > 0 || sys.S_Sol_Gen_ub_c[i] > 0
            Pa = res.P_sol_a[i] / 1000.0
            Pb = res.P_sol_b[i] / 1000.0
            Pc = res.P_sol_c[i] / 1000.0
            @printf("  Node %2d:  Pa=%.2f  Pb=%.2f  Pc=%.2f\n", i, Pa, Pb, Pc)
        end
    end
    
    if !isempty(vuf_buses)
        println("\n── Voltage Unbalance Factor (VUF) [%] ─────────────────────────")
        vuf = compute_vuf(res, sys.nN)
        for b in vuf_buses
            @printf("  Node %2d:  VUF = %.2f%%\n", b, vuf[b])
        end
    end
    
    println("\n" * "="^70 * "\n")
end

# ============================================================
# MAIN WORKFLOWS
# ============================================================

"""
    main1(; vuf_buses, vuf_threshold, n_oe_angles, save_plot, plot_file)

Workflow 1: Baseline OPF + Operating Envelope with fairness allocation.

Parameters:
- vuf_buses      : Nodes where VUF is constrained (empty = no constraint)
- vuf_threshold  : Maximum allowed VUF (default 0.02 = 2%)
- n_oe_angles    : Number of directional scans for OE boundary (default 36)
- save_plot      : Save figures to disk (default true)
- plot_file      : Output filename for OE plot
"""
function main1(;
    vuf_buses::Vector{Int}  = Int[],
    vuf_threshold::Float64  = 0.02,
    n_oe_angles::Int        = 36,
    save_plot::Bool         = true,
    plot_file::String       = "pq_operating_envelopes.png"
)
    println("\n" * "="^70)
    println("NETWORK 1 — WORKFLOW 1: OPF + Operating Envelopes")
    println("="^70)
    
    # Step 1: Build network
    println("\n[ 1/5 ]  Building system data …")
    sys = build_system_data()
    Y   = build_admittance_matrix(sys)

    # Step 2: Configure VUF constraints
    vuf_cfg = VUFConfig(vuf_buses; threshold = vuf_threshold)
    if isempty(vuf_buses)
        println("[ 2/5 ]  No VUF constraints requested.")
    else
        println("[ 2/5 ]  VUF ≤ $(vuf_threshold*100)% enforced at buses: $vuf_buses")
    end

    # Step 3: Solve baseline AC-OPF
    println("[ 3/5 ]  Solving baseline AC-OPF …")
    baseline_model = build_and_solve_opf(sys, Y, vuf_cfg)
    baseline_res   = extract_result(baseline_model, sys.nN)
    print_results(sys, baseline_res, vuf_buses)

    # Step 4: Compute operating envelopes
    println("[ 4/5 ]  Computing PQ Operating Envelopes …")
    oe_list = compute_operating_envelopes(sys, Y, vuf_cfg, baseline_res;
                                          n_angles = n_oe_angles)

    # Step 5: Plot and save
    println("[ 5/5 ]  Generating plots …")
    fig = plot_operating_envelopes(oe_list)
    if save_plot
        savefig(fig, plot_file)
        println("  Plot saved → $plot_file")
    end
    display(fig)

    return sys, baseline_res, oe_list, fig
end

"""
    main2(; vuf_buses, vuf_threshold, n_oe_angles, save_plot, plot_file)

Workflow 2: Baseline OPF + Dashboard visualization.

Generates four comprehensive figures:
  1. Active power (demand + PV)
  2. Reactive power (demand + PV)
  3. Voltage magnitudes + VUF
  4. Line flows + losses
"""
function main2(;
    vuf_buses::Vector{Int}  = Int[],
    vuf_threshold::Float64  = 0.02,
    n_oe_angles::Int        = 36,
    save_plot::Bool         = true,
    plot_dir::String        = "outputs/NW1_13Nodes/plots"
)
    println("\n" * "="^70)
    println("NETWORK 1 — WORKFLOW 2: OPF + Dashboard")
    println("="^70)
    
    # Step 1: Build network
    println("\n[ 1/4 ]  Building system data …")
    sys = build_system_data()
    Y   = build_admittance_matrix(sys)

    # Step 2: Configure VUF constraints
    vuf_cfg = VUFConfig(vuf_buses; threshold = vuf_threshold)
    if isempty(vuf_buses)
        println("[ 2/4 ]  No VUF constraints requested.")
    else
        println("[ 2/4 ]  VUF ≤ $(vuf_threshold*100)% enforced at buses: $vuf_buses")
    end

    # Step 3: Solve baseline AC-OPF
    println("[ 3/4 ]  Solving baseline AC-OPF …")
    baseline_model = build_and_solve_opf(sys, Y, vuf_cfg)
    baseline_res   = extract_result(baseline_model, sys.nN)
    print_results(sys, baseline_res, vuf_buses)

    # Step 4: Generate dashboard plots
    println("[ 4/4 ]  Generating dashboard …")
    save_dir = save_plot ? plot_dir : ""
    figs = plot_dashboard(sys, baseline_res; vuf_threshold, save_dir)

    return sys, baseline_res, figs
end

# ============================================================
# ENTRY POINT
# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    println("Three-Phase OPF — Network 1 (13 nodes)")
    
    # Run warm-up (triggers JIT compilation)
    println("\n(Compiling...)")
    main1(n_oe_angles = 1)
    main2()
    
    # Run full analysis
    println("\nRunning full analysis...")
    sys, res, oe, fig = main1(n_oe_angles = 36, vuf_buses = [5, 8, 11])
end
