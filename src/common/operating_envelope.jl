# ============================================================
# OPERATING ENVELOPE (OE) COMPUTATION
# ============================================================
# Trace the feasible P-Q region for each PV unit under network constraints.
# Uses directional scanning with sensitivity-weighted fairness allocation.

using JuMP, Ipopt, Printf
import MathOptInterface as MOI

include("structs.jl")
include("network_utils.jl")
include("optimization.jl")

"""
    find_pv_units(sys) -> Vector{PVUnit}

Identify all PV units in the system (nodes/phases with S_Sol_Gen_ub > 0).

Returns a vector of PVUnit structs, one per (node, phase) pair with capacity > 0.
"""
function find_pv_units(sys::SystemData)::Vector{PVUnit}
    units = PVUnit[]
    for i in 1:sys.nN
        sys.S_Sol_Gen_ub_a[i] > 0 && push!(units, PVUnit(i, :a))
        sys.S_Sol_Gen_ub_b[i] > 0 && push!(units, PVUnit(i, :b))
        sys.S_Sol_Gen_ub_c[i] > 0 && push!(units, PVUnit(i, :c))
    end
    return units
end

"""
    build_oe_model(sys, Y, vuf_cfg; fix_P, fix_Q, unit) -> JuMP.Model

Construct an OPF model for operating envelope tracing.

All PV units except the target unit are fixed to their baseline values.
This isolates the feasible region available to the target unit.

Parameters:
- fix_P, fix_Q : Dict of (node, phase) → value for fixed units
- unit         : Target PVUnit (included for completeness, not used directly)
"""
function build_oe_model(
    sys::SystemData, Y::Matrix{ComplexF64}, vuf_cfg::VUFConfig;
    fix_P::Dict{Tuple{Int,Symbol},Float64},
    fix_Q::Dict{Tuple{Int,Symbol},Float64},
    unit::PVUnit
)::JuMP.Model

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 200_000)
    set_optimizer_attribute(model, "tol",      1e-6)
    set_optimizer_attribute(model, "print_level", 0)
    set_optimizer_attribute(model, "linear_solver", "mumps")
    set_optimizer_attribute(model, "warm_start_init_point", "no")

    nN = sys.nN
    V0 = 400.0 / sqrt(3)

    # Variables
    @variable(model, P_G_a[1:nN]);  @variable(model, Q_G_a[1:nN])
    @variable(model, P_G_b[1:nN]);  @variable(model, Q_G_b[1:nN])
    @variable(model, P_G_c[1:nN]);  @variable(model, Q_G_c[1:nN])

    @variable(model, P_sol_a[1:nN]);  @variable(model, Q_sol_a[1:nN])
    @variable(model, P_sol_b[1:nN]);  @variable(model, Q_sol_b[1:nN])
    @variable(model, P_sol_c[1:nN]);  @variable(model, Q_sol_c[1:nN])

    @variable(model, V_a_r[1:nN], start = V0)
    @variable(model, V_a_i[1:nN], start = 0.0)
    @variable(model, V_b_r[1:nN], start = V0*cos(-2π/3))
    @variable(model, V_b_i[1:nN], start = V0*sin(-2π/3))
    @variable(model, V_c_r[1:nN], start = V0*cos( 2π/3))
    @variable(model, V_c_i[1:nN], start = V0*sin( 2π/3))

    # Slack bus reference (node 1, balanced)
    fix(V_a_r[1], V0*cos(0);     force=true)
    fix(V_a_i[1], V0*sin(0);     force=true)
    fix(V_b_r[1], V0*cos(-2π/3); force=true)
    fix(V_b_i[1], V0*sin(-2π/3); force=true)
    fix(V_c_r[1], V0*cos( 2π/3); force=true)
    fix(V_c_i[1], V0*sin( 2π/3); force=true)

    # Generation and PV limits
    for i in 1:nN
        s = sys.Gen_Status[i]
        for (PG, QG) in ((P_G_a,Q_G_a),(P_G_b,Q_G_b),(P_G_c,Q_G_c))
            @constraint(model, sys.P_Gen_lb[i]*s <= PG[i] <= sys.P_Gen_ub[i]*s)
            @constraint(model, sys.Q_Gen_lb[i]*s <= QG[i] <= sys.Q_Gen_ub[i]*s)
        end
    end
    @constraint(model,[i in 1:nN], P_sol_a[i]^2+Q_sol_a[i]^2 <= sys.S_Sol_Gen_ub_a[i]^2)
    @constraint(model,[i in 1:nN], P_sol_b[i]^2+Q_sol_b[i]^2 <= sys.S_Sol_Gen_ub_b[i]^2)
    @constraint(model,[i in 1:nN], P_sol_c[i]^2+Q_sol_c[i]^2 <= sys.S_Sol_Gen_ub_c[i]^2)

    # Voltage magnitude limits
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_a_r[i]^2+V_a_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_b_r[i]^2+V_b_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_c_r[i]^2+V_c_i[i]^2 <= sys.V_ub[i]^2)

    # Fix all non-target PV units at their baseline
    for ((node, ph), val) in fix_P
        var = (ph==:a) ? P_sol_a : (ph==:b) ? P_sol_b : P_sol_c
        fix(var[node], val; force=true)
    end
    for ((node, ph), val) in fix_Q
        var = (ph==:a) ? Q_sol_a : (ph==:b) ? Q_sol_b : Q_sol_c
        fix(var[node], val; force=true)
    end

    # VUF and power balance constraints
    add_vuf_constraints!(model, vuf_cfg, V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i)
    add_power_balance_exact!(model, sys, Y,
        V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i,
        P_G_a, P_G_b, P_G_c, Q_G_a, Q_G_b, Q_G_c,
        P_sol_a, P_sol_b, P_sol_c, Q_sol_a, Q_sol_b, Q_sol_c)

    return model
end

"""
    set_directional_objective!(model, unit, θ)

Update the OPF objective to maximize apparent power in direction θ.

Objective: max (cos(θ)·P + sin(θ)·Q) for the target PV unit.
This traces the OE boundary point in direction θ.
"""
function set_directional_objective!(model::JuMP.Model, unit::PVUnit, θ::Float64)
    P = unit.phase==:a ? model[:P_sol_a][unit.node] :
        unit.phase==:b ? model[:P_sol_b][unit.node] : model[:P_sol_c][unit.node]
    Q = unit.phase==:a ? model[:Q_sol_a][unit.node] :
        unit.phase==:b ? model[:Q_sol_b][unit.node] : model[:Q_sol_c][unit.node]
    @objective(model, Max, cos(θ)*P + sin(θ)*Q)
end

"""
    compute_operating_envelopes(sys, Y, vuf_cfg, baseline_res; n_angles) -> Vector{OEBoundary}

Compute the feasible P-Q region (operating envelope) for all PV units.

Algorithm:
  1. Identify all PV units (nodes/phases with capacity > 0)
  2. For each phase, fix all other units at their baseline values
  3. Perform directional scans at n_angles directions (0 to 2π)
  4. Use sensitivity-weighted fairness: allocate capacity proportional to
     each unit's voltage sensitivity
  5. Return OEBoundary for each unit

The resulting envelope respects all constraints:
  - AC power flow
  - Voltage limits
  - VUF constraints (if enabled)
  - PV capacity limits
"""
function compute_operating_envelopes(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int = 36
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units. Tracing OE ($n_angles angles, sensitivity-weighted fairness)...")

    # Helper functions for accessing baseline values
    get_P_base(u) = (u.phase==:a) ? baseline_res.P_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.P_sol_b[u.node] : baseline_res.P_sol_c[u.node]
    get_Q_base(u) = (u.phase==:a) ? baseline_res.Q_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.Q_sol_b[u.node] : baseline_res.Q_sol_c[u.node]

    # Pre-compute impedance matrix once
    Z_full = build_impedance_matrix(Y; slack_node=1)

    # Compute sensitivity weights for fair capacity allocation
    weights = compute_sensitivity_weights_cached(sys, baseline_res, units, Z_full)

    println("  Sensitivity weights (phase-normalised):")
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue
        println("    Phase $(uppercase(string(ph))):")
        for u in ph_units
            @printf("      Node %2d: w = %.4f\n", u.node, weights[u])
        end
    end

    # Angles for boundary scanning [0, 2π)
    angles  = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list = OEBoundary[]

    # Process each phase independently
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue

        println("  Phase $(uppercase(string(ph))) — $(length(ph_units)) units")

        # Fix all non-target phase units
        fix_P = Dict{Tuple{Int,Symbol},Float64}()
        fix_Q = Dict{Tuple{Int,Symbol},Float64}()
        for other in units
            if other.phase != ph
                fix_P[(other.node, other.phase)] = get_P_base(other)
                fix_Q[(other.node, other.phase)] = get_Q_base(other)
            end
        end

        # Build reusable OPF model for this phase
        model = build_oe_model(sys, Y, vuf_cfg;
                               fix_P=fix_P, fix_Q=fix_Q, unit=ph_units[1])

        # Slack variable α for sensitivity-weighted allocation
        @variable(model, alpha >= 0)
        @objective(model, Max, alpha)

        # Fairness constraints: each unit i expands by w_i·α in the given direction
        fair_p = Dict{Int, ConstraintRef}()
        fair_q = Dict{Int, ConstraintRef}()
        for u in ph_units
            Pvar = (u.phase==:a) ? model[:P_sol_a] :
                   (u.phase==:b) ? model[:P_sol_b] : model[:P_sol_c]
            Qvar = (u.phase==:a) ? model[:Q_sol_a] :
                   (u.phase==:b) ? model[:Q_sol_b] : model[:Q_sol_c]
            # Initial setup (θ=0): P = baseline + w·α, Q = baseline
            fair_p[u.node] = @constraint(model,
                Pvar[u.node] - weights[u]*alpha == get_P_base(u))
            fair_q[u.node] = @constraint(model,
                Qvar[u.node]                   == get_Q_base(u))
        end

        # Store boundary points for each unit
        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        # Directional scan
        for θ in angles
            dx = cos(θ)
            dy = sin(θ)

            # Update fairness coefficients for this direction
            for u in ph_units
                set_normalized_coefficient(fair_p[u.node], alpha, -weights[u] * dx)
                set_normalized_coefficient(fair_q[u.node], alpha, -weights[u] * dy)
            end

            # Solve for maximum α in this direction
            optimize!(model)

            ts = termination_status(model)
            if ts == MOI.LOCALLY_SOLVED || ts == MOI.OPTIMAL
                α_val = value(alpha)
                # Record boundary point: baseline + w·α·(cos θ, sin θ)
                for u in ph_units
                    push!(ph_bounds[u].P, get_P_base(u) + weights[u] * α_val * dx)
                    push!(ph_bounds[u].Q, get_Q_base(u) + weights[u] * α_val * dy)
                end
            else
                @printf("    ⚠  θ=%.1f°: %s\n", rad2deg(θ), ts)
            end
        end

        # Store results
        for u in ph_units
            n_pts = length(ph_bounds[u].P)
            @printf("    Node %2d (%s): %d/%d points, w=%.3f\n",
                    u.node, u.phase, n_pts, n_angles, weights[u])
            push!(oe_list, OEBoundary(u, ph_bounds[u].P, ph_bounds[u].Q,
                                      get_P_base(u), get_Q_base(u)))
        end
    end

    return oe_list
end

"""
    compute_sensitivity_weights_cached(sys, baseline_res, units, Z) -> Dict{PVUnit, Float64}

Compute sensitivity-based weights for PV capacity allocation.

Weight w_i for unit i = 1 / √(dV_i/dS), where dV_i/dS is the voltage
sensitivity (voltage rise per unit of apparent power injection).

Units with high sensitivity (weak bus) get lower weights. Units with low
sensitivity (strong bus) get higher weights, enabling fairer envelope tracing.

Accepts pre-computed impedance matrix Z to avoid redundant inversion.
"""
function compute_sensitivity_weights_cached(
    sys::SystemData,
    baseline_res::OPFResult,
    units::Vector{PVUnit},
    Z::Matrix{ComplexF64}
)::Dict{PVUnit, Float64}

    nN  = sys.nN
    raw = Dict{PVUnit, Float64}()

    # Build full 3N-element baseline voltage vector
    # Element 3*(i-1)+k corresponds to phase k of node i
    V_base_vec = Vector{ComplexF64}(undef, 3*nN)
    for i in 1:nN
        V_base_vec[3*(i-1)+1] = baseline_res.V_r_a[i] + im*baseline_res.V_i_a[i]
        V_base_vec[3*(i-1)+2] = baseline_res.V_r_b[i] + im*baseline_res.V_i_b[i]
        V_base_vec[3*(i-1)+3] = baseline_res.V_r_c[i] + im*baseline_res.V_i_c[i]
    end

    # Compute raw sensitivity for each unit
    for u in units
        phase_idx = (u.phase == :a) ? 1 : (u.phase == :b) ? 2 : 3
        col_idx = 3 * (u.node - 1) + phase_idx

        V_node = V_base_vec[col_idx]
        z_col = Z[:, col_idx]

        # Normalized impedance: Z / |V|
        z_norm = z_col ./ conj(V_node)

        # Sensitivity: √(dP² + dQ²) = √(Re² + Im²) of normalized impedance
        dVdP = norm(real.(z_norm))
        dVdQ = norm(imag.(z_norm))
        raw[u] = sqrt(dVdP^2 + dVdQ^2)
    end

    # Normalize within each phase for fairness
    result = Dict{PVUnit, Float64}()
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue

        raw_ph = [raw[u] for u in ph_units]
        max_raw = maximum(raw_ph)

        for u in ph_units
            result[u] = 1.0 / (raw[u] + 1e-10) / max(max_raw, 1e-10)
        end
    end

    return result
end
