# ============================================================
# OPERATING ENVELOPES
#   M1: sensitivity-weighted fairness   (compute_operating_envelopes)
#   M2: lexicographic max-min           (compute_operating_envelopes_lexmaxmin)
# ============================================================

function find_pv_units(sys::SystemData)::Vector{PVUnit}
    units = PVUnit[]
    for i in 1:sys.nN
        sys.S_Sol_Gen_ub_a[i] > 0 && push!(units, PVUnit(i, :a, sys.S_Sol_Gen_ub_a[i]))
        sys.S_Sol_Gen_ub_b[i] > 0 && push!(units, PVUnit(i, :b, sys.S_Sol_Gen_ub_b[i]))
        sys.S_Sol_Gen_ub_c[i] > 0 && push!(units, PVUnit(i, :c, sys.S_Sol_Gen_ub_c[i]))
    end
    return units
end

"""
    baseline_fixes(units, ph, baseline_res) -> (fix_P, fix_Q)

Fix every PV unit that is NOT on phase `ph` at its baseline injection.
"""
function baseline_fixes(units::Vector{PVUnit}, ph::Symbol, baseline_res::OPFResult)
    fix_P = Dict{Tuple{Int,Symbol},Float64}()
    fix_Q = Dict{Tuple{Int,Symbol},Float64}()
    for other in units
        if other.phase != ph
            fix_P[(other.node, other.phase)] = P_sol_of(baseline_res, other)
            fix_Q[(other.node, other.phase)] = Q_sol_of(baseline_res, other)
        end
    end
    return fix_P, fix_Q
end

"""
    compute_sensitivity_weights_cached(sys, baseline_res, units, Z; voltage_ref=:bus)

Phase-normalised weights w_u = 1/‖∂V/∂S_u‖ (scaled so max w = 1 per phase),
using a pre-computed Z = inv(Y).

`voltage_ref` selects the linearisation voltage:
  :bus  — divide each row k of Z[:,u] by conj(V_k⁰)   (exact linearisation)
  :unit — divide every row by conj(V_u⁰), the voltage at the unit's own bus
"""
function compute_sensitivity_weights_cached(
    sys::SystemData,
    baseline_res::OPFResult,
    units::Vector{PVUnit},
    Z::Matrix{ComplexF64};
    voltage_ref::Symbol = :bus
)::Dict{PVUnit, Float64}

    @assert voltage_ref in (:bus, :unit)  "voltage_ref must be :bus or :unit"
    nN  = sys.nN
    raw = Dict{PVUnit, Float64}()

    # Baseline voltage vector (3N), node-interleaved: 3(i-1)+1 = phase a of node i, …
    V_base_vec = Vector{ComplexF64}(undef, 3*nN)
    for i in 1:nN
        V_base_vec[3*(i-1)+1] = baseline_res.V_r_a[i] + im*baseline_res.V_i_a[i]
        V_base_vec[3*(i-1)+2] = baseline_res.V_r_b[i] + im*baseline_res.V_i_b[i]
        V_base_vec[3*(i-1)+3] = baseline_res.V_r_c[i] + im*baseline_res.V_i_c[i]
    end

    for u in units
        phase_idx = (u.phase == :a) ? 1 : (u.phase == :b) ? 2 : 3
        col_idx   = 3 * (u.node - 1) + phase_idx
        z_col     = Z[:, col_idx]

        V_ref = voltage_ref == :bus ? conj.(V_base_vec) : conj(V_base_vec[col_idx])
        dVdP  = real.(z_col ./ V_ref)
        dVdQ  = imag.(z_col ./ V_ref)

        σ      = norm(dVdP)^2 + norm(dVdQ)^2
        raw[u] = 1.0 / (sqrt(σ) + 1e-9)
    end

    # Phase-local normalisation
    weights = Dict{PVUnit, Float64}()
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue
        max_w = maximum(raw[u] for u in ph_units)
        for u in ph_units
            weights[u] = raw[u] / max_w
        end
    end

    return weights
end

"""
    compute_operating_envelopes(sys, Y, vuf_cfg, baseline_res; n_angles, voltage_ref)

M1 — sensitivity-weighted fairness. For each phase and direction θ, every unit u
on that phase moves by w_u·α·(cos θ, sin θ) and α is maximised. One model per
phase is reused across angles (only the α coefficients change).
"""
function compute_operating_envelopes(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int       = 36,
    voltage_ref::Symbol = :bus
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units. Tracing OE ($n_angles angles, sensitivity-weighted fairness)...")

    get_P_base(u) = P_sol_of(baseline_res, u)
    get_Q_base(u) = Q_sol_of(baseline_res, u)

    Z_full  = build_impedance_matrix(Y; slack_node=1)
    weights = compute_sensitivity_weights_cached(sys, baseline_res, units, Z_full;
                                                 voltage_ref=voltage_ref)

    println("  Sensitivity weights (phase-normalised):")
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue
        println("    Phase $(uppercase(string(ph))):")
        for u in ph_units
            @printf("      Node %2d: w = %.4f\n", u.node, weights[u])
        end
    end

    angles  = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list = OEBoundary[]

    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue

        println("  Phase $(uppercase(string(ph))) — $(length(ph_units)) units")

        fix_P, fix_Q = baseline_fixes(units, ph, baseline_res)
        model = build_base_model(sys, Y, vuf_cfg; fix_P=fix_P, fix_Q=fix_Q)

        @variable(model, alpha >= 0)
        @objective(model, Max, alpha)

        # Fairness constraints: unit i steps by w_i * alpha in the shared direction (cos θ, sin θ).
        # Coefficients on alpha are updated each angle via set_normalized_coefficient.
        fair_p = Dict{Int, ConstraintRef}()
        fair_q = Dict{Int, ConstraintRef}()
        for u in ph_units
            Pvar = _pick(u.phase, model[:P_sol_a], model[:P_sol_b], model[:P_sol_c])
            Qvar = _pick(u.phase, model[:Q_sol_a], model[:Q_sol_b], model[:Q_sol_c])
            # Initial coefficients (θ=0): P = P_base + w*alpha, Q = Q_base
            fair_p[u.node] = @constraint(model, Pvar[u.node] - weights[u]*alpha == get_P_base(u))
            fair_q[u.node] = @constraint(model, Qvar[u.node]                   == get_Q_base(u))
        end

        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        for θ in angles
            dx = cos(θ)
            dy = sin(θ)

            # Update only the alpha coefficient; RHS stays fixed at the baseline
            for u in ph_units
                set_normalized_coefficient(fair_p[u.node], alpha, -weights[u] * dx)
                set_normalized_coefficient(fair_q[u.node], alpha, -weights[u] * dy)
            end

            optimize!(model)

            ts = termination_status(model)
            if ts == MOI.LOCALLY_SOLVED || ts == MOI.OPTIMAL
                α_val = value(alpha)
                for u in ph_units
                    push!(ph_bounds[u].P, get_P_base(u) + weights[u] * α_val * dx)
                    push!(ph_bounds[u].Q, get_Q_base(u) + weights[u] * α_val * dy)
                end
            else
                # Record the baseline point (zero step), as M2 does, so every
                # unit keeps exactly one point per angle (needed for PoF alignment)
                @printf("    ⚠  θ=%.1f°: %s — recording baseline point\n", rad2deg(θ), ts)
                for u in ph_units
                    push!(ph_bounds[u].P, get_P_base(u))
                    push!(ph_bounds[u].Q, get_Q_base(u))
                end
            end
        end

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
    compute_operating_envelopes_lexmaxmin(sys, Y, vuf_cfg, baseline_res;
                                          n_angles, tol, stage_max_iter, final_max_iter)

M2 — pure lexicographic max-min. No sensitivity weights: the network alone
determines who binds. At each angle θ and stage:
- all free units move by the same scalar α in direction (cos θ, sin θ);
- α is maximised until the first network constraint binds;
- the unit(s) with the smallest apparent-power slack (within `tol`) are bound
  at that α and fixed; the rest proceed to the next stage.
If α stops growing (Δα < 1e-3) all remaining free units are bound at once.

`stage_max_iter` is the Ipopt iteration cap for intermediate stages (tol 1e-4);
`final_max_iter` applies to the last stage (tol 1e-6).
"""
function compute_operating_envelopes_lexmaxmin(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int       = 36,
    tol::Float64        = 1e-4,
    stage_max_iter::Int = 500,
    final_max_iter::Int = 2000
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units.")
    println("  Method: Pure Lexicographic Max-Min ($n_angles angles, tol=$tol)")

    get_P_base(u) = P_sol_of(baseline_res, u)
    get_Q_base(u) = Q_sol_of(baseline_res, u)

    angles  = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list = OEBoundary[]

    V_syms = (:V_a_r, :V_a_i, :V_b_r, :V_b_i, :V_c_r, :V_c_i)

    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue

        println("  Phase $(uppercase(string(ph))) — $(length(ph_units)) units")

        fix_P_base, fix_Q_base = baseline_fixes(units, ph, baseline_res)
        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        for θ in angles
            dx = cos(θ)
            dy = sin(θ)

            allocated  = Dict{PVUnit, Float64}()     # α at which each unit was bound
            free_units = copy(ph_units)
            warm       = Dict{Symbol, Vector{Float64}}()   # voltage warm start from previous stage

            stage  = 1
            α_prev = -Inf
            while !isempty(free_units)

                # Other phases at baseline + already-bound units at their achieved injection
                stage_fix_P = copy(fix_P_base)
                stage_fix_Q = copy(fix_Q_base)
                for (u, α_fixed) in allocated
                    stage_fix_P[(u.node, u.phase)] = get_P_base(u) + α_fixed * dx
                    stage_fix_Q[(u.node, u.phase)] = get_Q_base(u) + α_fixed * dy
                end

                model = build_base_model(sys, Y, vuf_cfg; fix_P=stage_fix_P, fix_Q=stage_fix_Q)

                # Relaxed tolerance for intermediate stages; full precision for the last one
                if length(free_units) > 1
                    set_optimizer_attribute(model, "tol",      1e-4)
                    set_optimizer_attribute(model, "max_iter", stage_max_iter)
                else
                    set_optimizer_attribute(model, "tol",      1e-6)
                    set_optimizer_attribute(model, "max_iter", final_max_iter)
                end

                if !isempty(warm)
                    for s in V_syms
                        set_start_value.(model[s], warm[s])
                    end
                    set_optimizer_attribute(model, "warm_start_init_point", "yes")
                end

                # Shared stage scalar — NO weight scaling
                @variable(model, α_stage >= 0)
                @objective(model, Max, α_stage)

                # Every free unit moves by exactly α in direction (dx, dy)
                for u in free_units
                    Pvar = _pick(u.phase, model[:P_sol_a], model[:P_sol_b], model[:P_sol_c])
                    Qvar = _pick(u.phase, model[:Q_sol_a], model[:Q_sol_b], model[:Q_sol_c])
                    @constraint(model, Pvar[u.node] - α_stage * dx == get_P_base(u))
                    @constraint(model, Qvar[u.node] - α_stage * dy == get_Q_base(u))
                end

                optimize!(model)

                ts = termination_status(model)
                if !(ts == MOI.LOCALLY_SOLVED || ts == MOI.OPTIMAL)
                    @printf("    ⚠  θ=%5.1f° stage %d: %s — zeroing remaining\n",
                            rad2deg(θ), stage, ts)
                    for u in free_units
                        allocated[u] = 0.0
                    end
                    break
                end

                for s in V_syms
                    warm[s] = value.(model[s])
                end

                α_val = value(α_stage)

                # α no longer grows → bind everyone still free
                if α_val - α_prev < 1e-3
                    @printf("    θ=%5.1f° stage %d: Δα=%.6f < tol — binding all free: %s\n",
                            rad2deg(θ), stage, α_val - α_prev,
                            join(["$(u.node)($(u.phase))" for u in free_units], ","))
                    for u in free_units
                        allocated[u] = α_val
                    end
                    break
                end
                α_prev = α_val

                # Bind the unit(s) with the smallest apparent-power slack at this α
                slacks = Dict{PVUnit, Float64}()
                for u in free_units
                    P_achieved = get_P_base(u) + α_val * dx
                    Q_achieved = get_Q_base(u) + α_val * dy
                    slacks[u]  = S_cap_of(sys, u)^2 - (P_achieved^2 + Q_achieved^2)
                end

                min_slack   = minimum(values(slacks))
                newly_bound = filter(u -> slacks[u] <= min_slack + tol, free_units)

                # Guard: always bind at least one unit per stage
                if isempty(newly_bound)
                    newly_bound = [free_units[argmin([slacks[u] for u in free_units])]]
                end

                for u in newly_bound
                    allocated[u] = α_val
                end
                filter!(u -> u ∉ newly_bound, free_units)

                @printf("    θ=%5.1f° stage %d: α*=%.5f | bound: %s | free: %s\n",
                        rad2deg(θ), stage, α_val,
                        join(["$(u.node)($(u.phase))" for u in newly_bound], ","),
                        join(["$(u.node)($(u.phase))" for u in free_units], ","))

                stage += 1
            end

            # Record boundary points — step is exactly α, no w_u
            for u in ph_units
                α_u = get(allocated, u, 0.0)
                push!(ph_bounds[u].P, get_P_base(u) + α_u * dx)
                push!(ph_bounds[u].Q, get_Q_base(u) + α_u * dy)
            end
        end

        for u in ph_units
            n_pts = length(ph_bounds[u].P)
            @printf("    Node %2d (%s): %d/%d boundary points\n",
                    u.node, u.phase, n_pts, n_angles)
            push!(oe_list, OEBoundary(u, ph_bounds[u].P, ph_bounds[u].Q,
                                      get_P_base(u), get_Q_base(u)))
        end
    end

    return oe_list
end
