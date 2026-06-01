# ============================================================
# NETWORK 2: IEEE 116-NODE DISTRIBUTION NETWORK
# ============================================================
# Three-phase unbalanced AC-OPF with dynamic PQ operating envelopes
# and voltage unbalance factor (VUF) constraints.
#
# Note: Full network data (topology, impedances, loads) extracted from
# the IEEE 123-node benchmark network and adapted to 116 nodes.

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
    get_z_matrix(r1, x1, r0, x0) -> Matrix{ComplexF64}

Build 3×3 phase impedance matrix from sequence parameters.

Standard utility power system convention:
  z1 = positive-sequence impedance [Ω]
  z0 = zero-sequence impedance [Ω]

Converts to phase components (abc) by eigenvalue decomposition.

Parameters:
- r1, x1 : Positive-sequence resistance and reactance [Ω/km]
- r0, x0 : Zero-sequence resistance and reactance [Ω/km]

Returns: 3×3 symmetric complex impedance matrix [Ω/km]
"""
function get_z_matrix(r1, x1, r0, x0)
    z1 = r1 + im*x1
    z0 = r0 + im*x0
    zs = (z0 + 2*z1) / 3
    zm = (z0 - z1) / 3
    Z = ComplexF64[
        zs  zm  zm;
        zm  zs  zm;
        zm  zm  zs]
    return Z
end

"""
    build_system_data() -> SystemData

Build the IEEE 116-node distribution network (derived from IEEE 123-node benchmark).

Network features:
- 116 nodes, 115 lines
- Longer lines with more realistic impedances (overhead + cable mix)
- More distributed unbalanced loads
- PV generators at 9 locations (3 per phase)
- Voltage limits: ±6% of nominal (225–264 V LN)

Returns a SystemData struct with all network parameters.
"""
function build_system_data()::SystemData
    nN = 116;  nL = 115
    
    # Initialize demand arrays
    P_Demand_a = zeros(Float64, nN); P_Demand_b = zeros(Float64, nN); P_Demand_c = zeros(Float64, nN)
    Q_Demand_a = zeros(Float64, nN); Q_Demand_b = zeros(Float64, nN); Q_Demand_c = zeros(Float64, nN)
    S_Sol_Gen_ub_a = zeros(Float64, nN); S_Sol_Gen_ub_b = zeros(Float64, nN); S_Sol_Gen_ub_c = zeros(Float64, nN)

    # Unbalanced demand pattern: Phase A heavily loaded, Phase C lightly loaded
    # Distributed across multiple nodes for realistic network behavior
    P_Demand_a[[5,11,21,29,34,46,50,57,64,72,76,81,88,103,107,113,116]] .= 3.0
    P_Demand_b[[7,12,25,30,37,42,47,52,66,73,77,83,93,97,104,110,114]] .= 2.0
    P_Demand_c[[10,13,26,33,39,45,49,54,62,67,75,87,95,102,106,112,115]] .= 1.0

    # Balanced 3-phase loads at selected nodes (e.g., large industrial/commercial sites)
    for bus in [40, 60, 80, 96]
        P_Demand_a[bus] = 2.5; P_Demand_b[bus] = 2.5; P_Demand_c[bus] = 2.5
        Q_Demand_a[bus] = 1.0;  Q_Demand_b[bus] = 1.0;  Q_Demand_c[bus] = 1.0
    end

    # Scale demand to realistic power levels [W]
    Scaling_factor = 1.25 * 1000
    P_Demand_a = P_Demand_a .* Scaling_factor
    P_Demand_b = P_Demand_b .* Scaling_factor
    P_Demand_c = P_Demand_c .* Scaling_factor
    Q_Demand_a = P_Demand_a .* 0.25
    Q_Demand_b = P_Demand_b .* 0.25
    Q_Demand_c = P_Demand_c .* 0.25

    # PV capacity per phase [VA] — distributed across network
    for b in [21, 34, 46, 57, 72, 80, 81, 103, 113]; S_Sol_Gen_ub_a[b] = 7500.0; end
    for b in [25, 37, 47, 66, 77, 80, 93, 104, 114]; S_Sol_Gen_ub_b[b] = 7500.0; end
    for b in [26, 39, 49, 62, 75, 80, 95, 106, 115]; S_Sol_Gen_ub_c[b] = 7500.0; end

    # Generation parameters (slack bus at node 1)
    P_Cost_G   = [1.0; zeros(nN-1)]
    P_Gen_lb   = [-5_000_000.0; zeros(nN-1)]
    P_Gen_ub   = [ 5_000_000.0; zeros(nN-1)]
    Q_Gen_lb   = [-1_000_000.0; zeros(nN-1)]
    Q_Gen_ub   = [ 1_000_000.0; zeros(nN-1)]
    Gen_Status = [1; zeros(Int,nN-1)]

    # Voltage limits [V]: ±6% of nominal (240 V LL → 139 V LN)
    V_base = 416.0 / sqrt(3)  # 240 V line-to-line
    V_lb   = fill(0.94 * V_base, nN)  # 225.76 V
    V_ub   = fill(1.06 * V_base, nN)  # 254.20 V

    # Network topology: from-bus and to-bus for 115 lines
    fbus = [
        1,2,2,3,3,4,4,6,6,8,8,9,9,14,14,15,15,16,16,17,17,18,31,31,36,36,19,19,20,20,
        22,22,23,23,24,27,27,28,28,32,32,35,35,38,38,44,44,41,41,43,43,48,48,51,51,61,
        61,58,58,53,53,55,55,56,56,59,59,63,63,65,65,78,78,69,69,70,68,68,71,71,74,74,
        79,79,98,98,105,105,82,82,84,84,85,85,86,89,89,111,111,90,90,91,91,92,92,94,94,
        99,100,101,99,108,108,109,109
    ]

    tbus = [
        2,3,4,5,10,6,8,7,13,9,14,11,12,15,19,16,28,17,37,31,18,36,34,39,40,42,20,21,22,
        23,25,33,26,24,27,29,30,32,49,35,46,38,48,41,44,47,50,43,45,52,53,51,54,58,61,
        62,64,60,68,55,76,56,59,57,67,63,65,66,72,69,78,80,81,75,70,77,71,97,73,74,79,
        83,82,98,103,105,106,107,84,85,87,95,88,86,89,90,111,112,115,91,92,93,102,94,99,
        96,104,100,101,108,109,110,114,113,116
    ]

    # Line impedance matrices from sequence parameters [Ω/km]
    # Different impedance classes based on conductor type and configuration
    
    # Class 1: Main feeder (larger conductors, lower impedance)
    Z1 = get_z_matrix(0.0793, 0.0518, 0.3144, 0.1640)
    # Class 2: Secondary feeder (medium conductors)
    Z2 = get_z_matrix(0.1528, 0.0877, 0.3744, 0.1885)
    # Class 3: Service lines (small conductors, higher impedance)
    Z3 = get_z_matrix(0.3368, 0.1935, 0.5976, 0.3310)

    # Line lengths [km] — mix of overhead and underground
    lengths = [
        0.50, 0.80, 0.85, 0.92, 1.10, 0.60, 1.50, 0.50, 0.92, 0.65, 1.20, 0.45, 0.92, 0.50,
        0.85, 0.45, 1.90, 0.50, 2.50, 0.50, 0.65, 0.92, 0.65, 1.20, 0.85, 0.92, 0.85, 1.10,
        0.50, 0.92, 0.65, 1.45, 0.50, 0.92, 0.65, 1.30, 0.85, 0.65, 0.92, 0.75, 0.50, 0.92,
        0.85, 0.50, 1.25, 0.92, 0.45, 0.75, 0.92, 0.85, 0.65, 0.92, 0.45, 0.65, 0.92, 0.85,
        0.50, 0.92, 0.75, 0.92, 0.85, 0.92, 0.75, 1.10, 0.50, 0.92, 1.25, 0.85, 0.92, 0.65,
        0.50, 0.75, 0.92, 0.50, 1.85, 0.92, 0.75, 0.92, 0.50, 0.92, 0.65, 0.85, 1.25, 0.75,
        0.85, 0.92, 1.50, 1.25, 0.92, 0.65, 0.85, 0.75, 0.92, 0.45, 0.65, 0.92, 0.85, 0.75,
        0.92, 0.85, 0.92, 0.75, 0.85, 0.92, 1.30, 0.92, 0.75, 0.92, 0.50, 0.85
    ]

    # Assign impedance class by line index (mix of classes for realism)
    Z_class = [
        repeat([Z1], 5)...,          # Lines 1-5: main feeder (class 1)
        repeat([Z2], 30)...,         # Lines 6-35: secondary feeder (class 2)
        repeat([Z3], 80)...,         # Lines 36-115: service lines (class 3)
    ]

    # Per-branch impedance matrices: Z_line = length × Z_class
    Zij = [lengths[k] .* Z_class[k] for k in 1:nL]

    # Shunt admittance (capacitive coupling, typically small in distribution)
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

Print formatted summary of OPF solution.
"""
function print_results(sys::SystemData, res::OPFResult, vuf_buses::Vector{Int}=Int[])
    println("\n" * "="^70)
    println("OPF SOLUTION SUMMARY (Network 2: 116-node IEEE)")
    println("="^70)
    
    if res.status != MOI.OPTIMAL && res.status != MOI.LOCALLY_SOLVED
        println("⚠  Solver status: $(res.status)")
        return
    end
    
    println("✓ Status: $(res.status)")
    println("  Objective value: $(round(res.obj; digits=2))")
    
    println("\n── Voltage Magnitudes [V] ─────────────────────────────────────")
    _mag(r, i) = sqrt(r^2 + i^2)
    for i in [1, 10, 50, 100, 116]
        i > sys.nN && continue
        Va = _mag(res.V_r_a[i], res.V_i_a[i])
        Vb = _mag(res.V_r_b[i], res.V_i_b[i])
        Vc = _mag(res.V_r_c[i], res.V_i_c[i])
        @printf("  Node %3d:  Va=%.1f  Vb=%.1f  Vc=%.1f\n", i, Va, Vb, Vc)
    end
    
    println("\n── Solar Injection [kW] ──────────────────────────────────────")
    pv_count = 0
    for i in 1:sys.nN
        if sys.S_Sol_Gen_ub_a[i] > 0 || sys.S_Sol_Gen_ub_b[i] > 0 || sys.S_Sol_Gen_ub_c[i] > 0
            Pa = res.P_sol_a[i] / 1000.0
            Pb = res.P_sol_b[i] / 1000.0
            Pc = res.P_sol_c[i] / 1000.0
            @printf("  Node %3d:  Pa=%.2f  Pb=%.2f  Pc=%.2f\n", i, Pa, Pb, Pc)
            pv_count += 1
            pv_count >= 5 && break
        end
    end
    
    if !isempty(vuf_buses)
        println("\n── Voltage Unbalance Factor (VUF) [%] ─────────────────────────")
        vuf = compute_vuf(res, sys.nN)
        for b in vuf_buses[1:min(5, length(vuf_buses))]
            @printf("  Node %3d:  VUF = %.2f%%\n", b, vuf[b])
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
"""
function main1(;
    vuf_buses::Vector{Int}  = Int[],
    vuf_threshold::Float64  = 0.02,
    n_oe_angles::Int        = 36,
    save_plot::Bool         = true,
    plot_file::String       = "pq_operating_envelopes_nw2.png"
)
    println("\n" * "="^70)
    println("NETWORK 2 — WORKFLOW 1: OPF + Operating Envelopes")
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
        println("[ 2/5 ]  VUF ≤ $(vuf_threshold*100)% enforced at $(length(vuf_buses)) buses")
    end

    # Step 3: Solve baseline AC-OPF
    println("[ 3/5 ]  Solving baseline AC-OPF (116 nodes, 115 lines) …")
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
    main2(; vuf_buses, vuf_threshold, save_plot, plot_dir)

Workflow 2: Baseline OPF + Dashboard visualization.
"""
function main2(;
    vuf_buses::Vector{Int}  = Int[],
    vuf_threshold::Float64  = 0.02,
    save_plot::Bool         = true,
    plot_dir::String        = "outputs/NW2_EULV_IEEE/plots"
)
    println("\n" * "="^70)
    println("NETWORK 2 — WORKFLOW 2: OPF + Dashboard")
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
        println("[ 2/4 ]  VUF ≤ $(vuf_threshold*100)% enforced at $(length(vuf_buses)) buses")
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
    println("Three-Phase OPF — Network 2 (116 nodes IEEE-based)")
    
    # Run warm-up (triggers JIT compilation)
    println("\n(Compiling...)")
    main1(n_oe_angles = 1)
    main2()
    
    # Run full analysis
    println("\nRunning full analysis...")
    sys, res, oe, fig = main1(n_oe_angles = 36)
end
