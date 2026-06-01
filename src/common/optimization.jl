# ============================================================
# OPTIMIZATION AND CONSTRAINT FUNCTIONS
# ============================================================
# Core OPF model builder and solver, with constraint definitions for:
# - Nodal power balance (AC) via exact complex formulation
# - Voltage unbalance factor (VUF) limits
# - Voltage magnitude bounds
# - Power generation and solar injection limits

using JuMP, Ipopt
import MathOptInterface as MOI

include("structs.jl")

"""
    add_vuf_constraints!(model, vuf_cfg, V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i)

Add Voltage Unbalance Factor (VUF) constraints to the OPF model.

For each bus in vuf_cfg.buses, enforce:
  1. |V₊|² ≥ pos_seq_floor (prevents trivial solutions with |V+| → 0)
  2. |V₋|² ≤ ε² · |V₊|² where ε = vuf_cfg.threshold

Uses Clarke-Concordia decomposition (no auxiliary variables required).

Parameters:
- model          : JuMP optimization model
- vuf_cfg        : VUFConfig specifying buses and thresholds
- V_*_r, V_*_i   : Voltage real/imaginary components per node
"""
function add_vuf_constraints!(
    model::JuMP.Model,
    vuf_cfg::VUFConfig,
    V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i
)
    isempty(vuf_cfg.buses) && return
    
    # Clarke-Concordia matrix factors for 120° phase shift: a = e^(j·2π/3)
    ar, ai  = -0.5,  sqrt(3)/2
    a2r, a2i = -0.5, -sqrt(3)/2
    ε2  = vuf_cfg.threshold^2
    flr = vuf_cfg.pos_seq_floor

    for b in vuf_cfg.buses
        # Positive sequence: V+ = (1/3)·(Va + a·Vb + a²·Vc)
        Vp_r = @expression(model, (1/3)*( V_a_r[b] + ar*V_b_r[b] - ai*V_b_i[b] + a2r*V_c_r[b] - a2i*V_c_i[b] ))
        Vp_i = @expression(model, (1/3)*( V_a_i[b] + ar*V_b_i[b] + ai*V_b_r[b] + a2r*V_c_i[b] + a2i*V_c_r[b] ))
        
        # Negative sequence: V- = (1/3)·(Va + a²·Vb + a·Vc)
        Vn_r = @expression(model, (1/3)*( V_a_r[b] + a2r*V_b_r[b] - a2i*V_b_i[b] + ar*V_c_r[b] - ai*V_c_i[b] ))
        Vn_i = @expression(model, (1/3)*( V_a_i[b] + a2r*V_b_i[b] + a2i*V_b_r[b] + ar*V_c_i[b] + ai*V_c_r[b] ))

        # Floor on |V+|² and VUF inequality
        @constraint(model, Vp_r^2 + Vp_i^2 >= flr)
        @constraint(model, Vn_r^2 + Vn_i^2 <= ε2 * (Vp_r^2 + Vp_i^2))
    end
end

"""
    add_power_balance_exact!(model, sys, Y, V_*_r, V_*_i, ...)

Add nodal AC power balance constraints to the OPF model.

Uses the exact three-phase load flow formulation:
  S = diag(V) · (Y · conj(V)) ⟹ P + jQ per phase

Splits Y = G + jB (real and imaginary parts) and imposes:
  G[i] + j·Q[i] = (PG[i] + Psol[i] - Pdemand[i]) + j(QG[i] + Qsol[i] - Qdemand[i])

Parameters:
- model          : JuMP optimization model
- sys            : SystemData with network parameters
- Y              : Admittance matrix (built via build_admittance_matrix)
- V_*_r, V_*_i   : Voltage real/imaginary components per node
- P_G_*, Q_G_*   : Generation per phase
- P_sol_*, Q_sol_*: Solar injection per phase
"""
function add_power_balance_exact!(
    model::JuMP.Model,
    sys::SystemData,
    Y::Matrix{ComplexF64},
    V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i,
    P_G_a, P_G_b, P_G_c, Q_G_a, Q_G_b, Q_G_c,
    P_sol_a, P_sol_b, P_sol_c, Q_sol_a, Q_sol_b, Q_sol_c
)
    nN = sys.nN
    n3 = 3 * nN

    # Split admittance into real and imaginary parts
    G = real.(Y)
    B = imag.(Y)

    # Helper functions to map stacked index k ∈ [1, 3nN] to (phase, node)
    @inline function vr(k::Int)
        p = (k - 1) % 3 + 1
        i = (k - 1) ÷ 3 + 1
        return p==1 ? V_a_r[i] : p==2 ? V_b_r[i] : V_c_r[i]
    end
    @inline function vi(k::Int)
        p = (k - 1) % 3 + 1
        i = (k - 1) ÷ 3 + 1
        return p==1 ? V_a_i[i] : p==2 ? V_b_i[i] : V_c_i[i]
    end

    # Current injection: J = Y · conj(V)
    #   J_r[k] = Σ_j ( G[k,j]·V_r[j] + B[k,j]·V_i[j] )
    #   J_i[k] = Σ_j ( B[k,j]·V_r[j] - G[k,j]·V_i[j] )
    @expression(model, J_r[k=1:n3], sum(G[k,j]*vr(j) + B[k,j]*vi(j) for j in 1:n3))
    @expression(model, J_i[k=1:n3], sum(B[k,j]*vr(j) - G[k,j]*vi(j) for j in 1:n3))

    # Power: S = V ⊙ J (element-wise product with conjugate current already in J)
    #   P[k] = V_r[k]·J_r[k] - V_i[k]·J_i[k]
    #   Q[k] = V_r[k]·J_i[k] + V_i[k]·J_r[k]
    @expression(model, P[k=1:n3], vr(k)*J_r[k] - vi(k)*J_i[k])
    @expression(model, Q[k=1:n3], vr(k)*J_i[k] + vi(k)*J_r[k])

    # Map back to (node, phase) and impose power balance at each node
    for i in 1:nN
        ka, kb, kc = 3i-2, 3i-1, 3i

        # Active power balance
        @constraint(model, P_G_a[i] + P_sol_a[i] - sys.P_Demand_a[i] == P[ka])
        @constraint(model, P_G_b[i] + P_sol_b[i] - sys.P_Demand_b[i] == P[kb])
        @constraint(model, P_G_c[i] + P_sol_c[i] - sys.P_Demand_c[i] == P[kc])

        # Reactive power balance
        @constraint(model, Q_G_a[i] + Q_sol_a[i] - sys.Q_Demand_a[i] == Q[ka])
        @constraint(model, Q_G_b[i] + Q_sol_b[i] - sys.Q_Demand_b[i] == Q[kb])
        @constraint(model, Q_G_c[i] + Q_sol_c[i] - sys.Q_Demand_c[i] == Q[kc])
    end
end

"""
    build_and_solve_opf(sys, Y, vuf_cfg; ...) -> JuMP.Model

Build and solve the three-phase AC Optimal Power Flow (OPF) model.

Features:
- Exact AC power flow via rectangular voltage coordinates
- Optional VUF constraints (demand-side voltage regulation)
- Optional directional objective (for operating envelope scans)
- Optional fixed PV setpoints (for OE exploration)

Parameters:
- sys                 : SystemData
- Y                   : Admittance matrix (from build_admittance_matrix)
- vuf_cfg             : VUFConfig (constraints on unbalance)
- fix_P_sol, fix_Q_sol: Dict of fixed PV injections (node, phase) → value
- obj_r_node, obj_r_phase, obj_angle : Objective targeting (OE computation)
- maximise            : Direction (min if false, max if true)

Returns: Optimized JuMP.Model
"""
function build_and_solve_opf(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig            = VUFConfig();
    fix_P_sol::Dict{Tuple{Int,Symbol},Float64} = Dict{Tuple{Int,Symbol},Float64}(),
    fix_Q_sol::Dict{Tuple{Int,Symbol},Float64} = Dict{Tuple{Int,Symbol},Float64}(),
    obj_r_node::Int      = 0,
    obj_r_phase::Symbol  = :a,
    obj_angle::Float64   = 0.0,
    maximise::Bool       = false
)::JuMP.Model

    model = Model(Ipopt.Optimizer)
    
    # Conservative solver tolerances for robust convergence
    set_optimizer_attribute(model, "max_iter", 200_000)
    set_optimizer_attribute(model, "tol",      1e-6)
    set_optimizer_attribute(model, "print_level", 0)
    set_optimizer_attribute(model, "linear_solver", "mumps")
    set_optimizer_attribute(model, "warm_start_init_point", "no")

    nN = sys.nN
    V0 = 400.0 / sqrt(3)  # Nominal LN voltage (230V)

    # ──────────────────────────────────────────────────────────────────────
    # Decision Variables
    # ──────────────────────────────────────────────────────────────────────
    
    # Generation per phase
    @variable(model, P_G_a[1:nN]);  @variable(model, Q_G_a[1:nN])
    @variable(model, P_G_b[1:nN]);  @variable(model, Q_G_b[1:nN])
    @variable(model, P_G_c[1:nN]);  @variable(model, Q_G_c[1:nN])

    # Solar injection per phase
    @variable(model, P_sol_a[1:nN]);  @variable(model, Q_sol_a[1:nN])
    @variable(model, P_sol_b[1:nN]);  @variable(model, Q_sol_b[1:nN])
    @variable(model, P_sol_c[1:nN]);  @variable(model, Q_sol_c[1:nN])

    # Voltage (rectangular coordinates, per phase)
    @variable(model, V_a_r[1:nN], start = V0)
    @variable(model, V_a_i[1:nN], start = 0.0)
    @variable(model, V_b_r[1:nN], start = V0*cos(-2π/3))
    @variable(model, V_b_i[1:nN], start = V0*sin(-2π/3))
    @variable(model, V_c_r[1:nN], start = V0*cos( 2π/3))
    @variable(model, V_c_i[1:nN], start = V0*sin( 2π/3))

    # ──────────────────────────────────────────────────────────────────────
    # Objective
    # ──────────────────────────────────────────────────────────────────────
    
    if obj_r_node > 0
        # Directional scan: optimize one PV unit in a specified direction
        P_tgt = (obj_r_phase==:a) ? P_sol_a[obj_r_node] :
                (obj_r_phase==:b) ? P_sol_b[obj_r_node] : P_sol_c[obj_r_node]
        Q_tgt = (obj_r_phase==:a) ? Q_sol_a[obj_r_node] :
                (obj_r_phase==:b) ? Q_sol_b[obj_r_node] : Q_sol_c[obj_r_node]
        sense = maximise ? -1.0 : 1.0
        @objective(model, Min, sense*(cos(obj_angle)*P_tgt + sin(obj_angle)*Q_tgt))
    else
        # Default: minimize generation cost
        @objective(model, Min,
            sum(sys.P_Cost_G[i]*(P_G_a[i]+P_G_b[i]+P_G_c[i]) for i in 1:nN))
    end

    # ──────────────────────────────────────────────────────────────────────
    # Slack Bus (node 1) — Fixed balanced reference
    # ──────────────────────────────────────────────────────────────────────
    
    fix(V_a_r[1], V0*cos(0);     force=true)
    fix(V_a_i[1], V0*sin(0);     force=true)
    fix(V_b_r[1], V0*cos(-2π/3); force=true)
    fix(V_b_i[1], V0*sin(-2π/3); force=true)
    fix(V_c_r[1], V0*cos( 2π/3); force=true)
    fix(V_c_i[1], V0*sin( 2π/3); force=true)

    # ──────────────────────────────────────────────────────────────────────
    # Generation Limits
    # ──────────────────────────────────────────────────────────────────────
    
    for i in 1:nN
        s = sys.Gen_Status[i]
        for (PG, QG) in ((P_G_a,Q_G_a),(P_G_b,Q_G_b),(P_G_c,Q_G_c))
            @constraint(model, sys.P_Gen_lb[i]*s <= PG[i] <= sys.P_Gen_ub[i]*s)
            @constraint(model, sys.Q_Gen_lb[i]*s <= QG[i] <= sys.Q_Gen_ub[i]*s)
        end
    end

    # ──────────────────────────────────────────────────────────────────────
    # PV Capacity Limits (apparent power)
    # ──────────────────────────────────────────────────────────────────────
    
    @constraint(model,[i in 1:nN], P_sol_a[i]^2+Q_sol_a[i]^2 <= sys.S_Sol_Gen_ub_a[i]^2)
    @constraint(model,[i in 1:nN], P_sol_b[i]^2+Q_sol_b[i]^2 <= sys.S_Sol_Gen_ub_b[i]^2)
    @constraint(model,[i in 1:nN], P_sol_c[i]^2+Q_sol_c[i]^2 <= sys.S_Sol_Gen_ub_c[i]^2)

    # ──────────────────────────────────────────────────────────────────────
    # Fixed PV Setpoints (for envelope boundary scans)
    # ──────────────────────────────────────────────────────────────────────
    
    for ((node, ph), val) in fix_P_sol
        var = (ph==:a) ? P_sol_a : (ph==:b) ? P_sol_b : P_sol_c
        fix(var[node], val; force=true)
    end
    for ((node, ph), val) in fix_Q_sol
        var = (ph==:a) ? Q_sol_a : (ph==:b) ? Q_sol_b : Q_sol_c
        fix(var[node], val; force=true)
    end

    # ──────────────────────────────────────────────────────────────────────
    # Voltage Magnitude Limits
    # ──────────────────────────────────────────────────────────────────────
    
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_a_r[i]^2+V_a_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_b_r[i]^2+V_b_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_c_r[i]^2+V_c_i[i]^2 <= sys.V_ub[i]^2)

    # ──────────────────────────────────────────────────────────────────────
    # VUF Constraints (optional, per-bus)
    # ──────────────────────────────────────────────────────────────────────
    
    add_vuf_constraints!(model, vuf_cfg, V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i)

    # ──────────────────────────────────────────────────────────────────────
    # AC Power Balance Constraints (exact three-phase)
    # ──────────────────────────────────────────────────────────────────────
    
    add_power_balance_exact!(model, sys, Y,
        V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i,
        P_G_a, P_G_b, P_G_c, Q_G_a, Q_G_b, Q_G_c,
        P_sol_a, P_sol_b, P_sol_c, Q_sol_a, Q_sol_b, Q_sol_c)

    # Solve
    optimize!(model)
    return model
end

"""
    extract_result(model, nN) -> OPFResult

Extract optimized solution from the JuMP model into an OPFResult struct.

Captures all voltage and generation variables in a decoupled form,
suitable for downstream post-processing and visualization.
"""
function extract_result(model::JuMP.Model, nN::Int)::OPFResult
    ts = termination_status(model)
    obj = (ts == MOI.LOCALLY_SOLVED || ts == MOI.OPTIMAL) ? objective_value(model) : NaN
    return OPFResult(ts, obj,
        value.(model[:V_a_r]), value.(model[:V_a_i]),
        value.(model[:V_b_r]), value.(model[:V_b_i]),
        value.(model[:V_c_r]), value.(model[:V_c_i]),
        value.(model[:P_G_a]), value.(model[:Q_G_a]),
        value.(model[:P_G_b]), value.(model[:Q_G_b]),
        value.(model[:P_G_c]), value.(model[:Q_G_c]),
        value.(model[:P_sol_a]), value.(model[:Q_sol_a]),
        value.(model[:P_sol_b]), value.(model[:Q_sol_b]),
        value.(model[:P_sol_c]), value.(model[:Q_sol_c]))
end
