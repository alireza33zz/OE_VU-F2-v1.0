# ============================================================
# UTILITARIAN BENCHMARK + PRICE OF FAIRNESS
# ============================================================
#
# Same feasible set F, baseline and direction θ_m as M1/M2, but each unit gets
# its own step δ_u and the objective is max Σ_u δ_u. M1 (δ_u = w_u α*) and
# M2 (δ_u = α_u) are both feasible points of this problem, so
#
#   PoF^k = 1 − Σ_u A_u^k / Σ_u A_u^util      (A_u = envelope area of unit u)
#
# U* is a locally optimal benchmark (nonconvex NLP), warm-started from the
# better of the two fair solutions at each angle.
# ============================================================

# P/Q JuMP variables of a unit
_pq_vars(model::JuMP.Model, u::PVUnit) =
    (_pick(u.phase, model[:P_sol_a], model[:P_sol_b], model[:P_sol_c])[u.node],
     _pick(u.phase, model[:Q_sol_a], model[:Q_sol_b], model[:Q_sol_c])[u.node])

# Step length of an OE boundary point from the baseline
_step(oe::OEBoundary, m::Int) =
    hypot(oe.P_boundary[m] - oe.P_baseline, oe.Q_boundary[m] - oe.Q_baseline)

"""
    compute_operating_envelopes_utilitarian(sys, Y, vuf_cfg, baseline_res;
                                            n_angles, warm_lists) -> (oe_list, U_star)

Utilitarian benchmark: per-unit steps δ_u ≥ 0 along the shared direction θ_m,
objective max Σ_u δ_u. `warm_lists = [oe_M1, oe_M2]`: at each angle the fair
solution with the larger aggregate step is used as starting point.
`U_star[(phase, m)] = Σ_u δ_u*` (NaN where the solve failed).
"""
function compute_operating_envelopes_utilitarian(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int = 36,
    warm_lists::Vector{Vector{OEBoundary}} = Vector{OEBoundary}[]
)
    units = find_pv_units(sys)
    get_P_base(u) = P_sol_of(baseline_res, u)
    get_Q_base(u) = Q_sol_of(baseline_res, u)

    warm_maps = [Dict((oe.unit.node, oe.unit.phase) => oe for oe in L) for L in warm_lists]
    angles    = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list   = OEBoundary[]
    U_star    = Dict{Tuple{Symbol,Int},Float64}()

    println("  Utilitarian benchmark ($n_angles angles)")
    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue
        n_ph = length(ph_units)

        println("  Phase $(uppercase(string(ph))) — $n_ph units")

        fix_P, fix_Q = baseline_fixes(units, ph, baseline_res)
        model = build_base_model(sys, Y, vuf_cfg; fix_P=fix_P, fix_Q=fix_Q)

        @variable(model, δ[1:n_ph] >= 0)
        @objective(model, Max, sum(δ))

        # P_u = P̂_u + δ_u cosθ ,  Q_u = Q̂_u + δ_u sinθ  (coefficients updated per angle)
        cp = Vector{ConstraintRef}(undef, n_ph)
        cq = Vector{ConstraintRef}(undef, n_ph)
        for (j, u) in enumerate(ph_units)
            Pv, Qv = _pq_vars(model, u)
            cp[j] = @constraint(model, Pv - δ[j] == get_P_base(u))
            cq[j] = @constraint(model, Qv        == get_Q_base(u))
        end

        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        for (m, θ) in enumerate(angles)
            dx, dy = cos(θ), sin(θ)
            for j in 1:n_ph
                set_normalized_coefficient(cp[j], δ[j], -dx)
                set_normalized_coefficient(cq[j], δ[j], -dy)
            end

            # Warm start from the better fair solution at this angle
            best, bestU = nothing, -Inf
            for mp in warm_maps
                oes = [mp[(u.node, u.phase)] for u in ph_units]
                all(length(o.P_boundary) == n_angles for o in oes) || continue
                steps = [_step(o, m) for o in oes]
                if sum(steps) > bestU
                    best, bestU = steps, sum(steps)
                end
            end
            if best !== nothing
                for (j, u) in enumerate(ph_units)
                    Pv, Qv = _pq_vars(model, u)
                    set_start_value(δ[j], best[j])
                    set_start_value(Pv, get_P_base(u) + best[j]*dx)
                    set_start_value(Qv, get_Q_base(u) + best[j]*dy)
                end
            end

            optimize!(model)
            ts = termination_status(model)
            if ts == MOI.LOCALLY_SOLVED || ts == MOI.OPTIMAL
                δv = value.(δ)
                U_star[(ph, m)] = sum(δv)
                for (j, u) in enumerate(ph_units)
                    push!(ph_bounds[u].P, get_P_base(u) + δv[j]*dx)
                    push!(ph_bounds[u].Q, get_Q_base(u) + δv[j]*dy)
                end
            else
                @printf("    ⚠  θ=%.1f°: %s (benchmark point = baseline)\n", rad2deg(θ), ts)
                U_star[(ph, m)] = NaN
                for u in ph_units           # keep index alignment
                    push!(ph_bounds[u].P, get_P_base(u))
                    push!(ph_bounds[u].Q, get_Q_base(u))
                end
            end
        end

        for u in ph_units
            push!(oe_list, OEBoundary(u, ph_bounds[u].P, ph_bounds[u].Q,
                                      get_P_base(u), get_Q_base(u)))
        end
    end
    return oe_list, U_star
end

"""
    aggregate_steps(oe_list, n_angles) -> Dict((phase, m) => Σ_u step)
"""
function aggregate_steps(oe_list::Vector{OEBoundary}, n_angles::Int)
    U = Dict{Tuple{Symbol,Int},Float64}()
    for oe in oe_list
        length(oe.P_boundary) == n_angles ||
            error("Node $(oe.unit.node)($(oe.unit.phase)): $(length(oe.P_boundary))/$n_angles points — failed angles break alignment")
        for m in 1:n_angles
            U[(oe.unit.phase, m)] = get(U, (oe.unit.phase, m), 0.0) + _step(oe, m)
        end
    end
    return U
end

# Shoelace area (points are in angular trace order around the baseline)
function polygon_area(oe::OEBoundary)
    x, y = oe.P_boundary, oe.Q_boundary
    n = length(x); n < 3 && return 0.0
    return 0.5 * abs(sum(x[i]*y[mod1(i+1,n)] - x[mod1(i+1,n)]*y[i] for i in 1:n))
end

"""
    price_of_fairness(oe_fair, oe_util, U_star, n_angles) -> Float64

PoF = 1 − ΣA_fair / ΣA_util  (total envelope area of the fair method vs the
utilitarian benchmark, summed over all units).

Sanity check: warns if, at some (phase, angle), the benchmark's aggregate step
Σδ_u ended below the fair one — a local optimum of the utilitarian NLP.
"""
function price_of_fairness(oe_fair::Vector{OEBoundary}, oe_util::Vector{OEBoundary},
                           U_star::Dict{Tuple{Symbol,Int},Float64}, n_angles::Int)
    U_fair  = aggregate_steps(oe_fair, n_angles)
    n_local = count(k -> !isnan(U_star[k]) && U_fair[k] > U_star[k] + 1e-3, keys(U_star))
    n_local > 0 && @warn "Utilitarian benchmark below the fair solution at $n_local (phase, angle) pairs — local optimum"

    return 1 - sum(polygon_area, oe_fair) / sum(polygon_area, oe_util)
end
