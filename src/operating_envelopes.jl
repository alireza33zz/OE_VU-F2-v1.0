# ============================================================
# operating_envelopes.jl — OE computation (weighted fairness
#                           and lexicographic max-min)
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

angular_diff(b, a) = mod(b - a + π, 2π) - π

"""
    local_optimal_angle(u, θ_global, weights, Z, baseline_res, sys)

Rotate the global tracing angle θ_global toward the node's locally-optimal
P–Q direction, blended by its sensitivity weight.

High-sensitivity nodes (small w) are rotated more aggressively toward their
optimal reactive direction; low-sensitivity nodes (w ≈ 1) follow θ_global.
"""
function local_optimal_angle(
    u::PVUnit,
    θ_global::Float64,
    weights::Dict{PVUnit,Float64},
    Z::Matrix{ComplexF64},
    baseline_res::OPFResult,
    sys::SystemData
)::Float64

    nN = sys.nN
    phase_idx = (u.phase == :a) ? 1 : (u.phase == :b) ? 2 : 3
    col_idx = 3 * (u.node - 1) + phase_idx

    V_base = (u.phase == :a) ? (baseline_res.V_r_a[u.node] + im*baseline_res.V_i_a[u.node]) :
             (u.phase == :b) ? (baseline_res.V_r_b[u.node] + im*baseline_res.V_i_b[u.node]) :
                               (baseline_res.V_r_c[u.node] + im*baseline_res.V_i_c[u.node])

    z_col = Z[:, col_idx]
    dVdP  = norm(real.(z_col ./ conj(V_base)))
    dVdQ  = norm(imag.(z_col ./ conj(V_base)))

    θ_opt = atan(dVdQ, dVdP) + π
    blend = 1.0 - weights[u]
    θ_i   = θ_global + blend * angular_diff(θ_opt, θ_global)

    return θ_i
end

"""
    compute_sensitivity_weights_cached(sys, baseline_res, units, Z)

Same as compute_sensitivity_weights but accepts a pre-computed Z = inv(Y).
Call this version inside compute_operating_envelopes to avoid re-inverting Y.
"""
function compute_sensitivity_weights_cached(
    sys::SystemData,
    baseline_res::OPFResult,
    units::Vector{PVUnit},
    Z::Matrix{ComplexF64}
)::Dict{PVUnit, Float64}

    nN  = sys.nN
    raw = Dict{PVUnit, Float64}()

    V_base_vec = Vector{ComplexF64}(undef, 3*nN)
    for i in 1:nN
        V_base_vec[3*(i-1)+1] = baseline_res.V_r_a[i] + im*baseline_res.V_i_a[i]
        V_base_vec[3*(i-1)+2] = baseline_res.V_r_b[i] + im*baseline_res.V_i_b[i]
        V_base_vec[3*(i-1)+3] = baseline_res.V_r_c[i] + im*baseline_res.V_i_c[i]
    end

    for u in units
        phase_idx = (u.phase == :a) ? 1 : (u.phase == :b) ? 2 : 3
        col_idx   = 3 * (u.node - 1) + phase_idx

        z_col = Z[:, col_idx]

        # Correct: divide each element k by conj(V_k^0), not by conj(V_u^0)
        dVdP = real.(z_col ./ conj.(V_base_vec))
        dVdQ = imag.(z_col ./ conj.(V_base_vec))

        σ      = norm(dVdP)^2 + norm(dVdQ)^2
        raw[u] = 1.0 / (sqrt(σ) + 1e-9)
    end

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
    plot_operating_envelopes(oe_list) -> Plots.Plot
"""
function plot_operating_envelopes(oe_list::Vector{OEBoundary})
    n     = length(oe_list)
    nrows = floor(Int, sqrt(n * 1.0))
    ncols = ceil(Int, n / nrows)
    scale = 1e-3   # W → kW

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=900)

    θ_circ = range(0, 2π, length=200)

    for (k, oe) in enumerate(oe_list)
        ax    = p[k]
        S_max = oe.unit.S_max * scale

        plot!(ax, S_max.*cos.(θ_circ), S_max.*sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label="S_max")

        if length(oe.P_boundary) >= 3
            cx = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
            idx = sortperm(atan.(oe.Q_boundary.-cy, oe.P_boundary.-cx))
            Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
            Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale
            plot!(ax, Ps, Qs; seriestype=:shape,
                  fillcolor=:dodgerblue, fillalpha=0.18,
                  lc=:dodgerblue, lw=1.8, label="OE")
        end

        scatter!(ax, [oe.P_baseline*scale], [oe.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label="Baseline")

        hline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)
        vline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)

        title!(ax, "Node $(oe.unit.node) · Ph $(uppercase(string(oe.unit.phase)))")
        xlabel!(ax, "P (kW)");  ylabel!(ax, "Q (kvar)")
        plot!(ax, legend=(k==1), legendfont=font(7), titlefont=font(9,:bold),
              tickfont=font(7), guidefont=font(8),
              aspect_ratio=:equal, framestyle=:box, grid=true, gridalpha=0.25)
    end

    # Blank unused grid cells
    for k in n+1:nrows*ncols
        plot!(p[k], axis=false, ticks=false, grid=false, label=false)
    end
    return p
end

# ── Weighted-fairness OE ────────────────────────────────────
"""
    compute_operating_envelopes(sys, Y, vuf_cfg, baseline_res; n_angles, V_base)
"""
function compute_operating_envelopes(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int    = 36
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units. Tracing OE ($n_angles angles, sensitivity-weighted fairness)...")

    get_P_base(u) = (u.phase==:a) ? baseline_res.P_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.P_sol_b[u.node] : baseline_res.P_sol_c[u.node]
    get_Q_base(u) = (u.phase==:a) ? baseline_res.Q_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.Q_sol_b[u.node] : baseline_res.Q_sol_c[u.node]

    Z_full  = build_impedance_matrix(Y; slack_node=1)
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

    angles  = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list = OEBoundary[]

    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue

        println("  Phase $(uppercase(string(ph))) — $(length(ph_units)) units")

        fix_P = Dict{Tuple{Int,Symbol},Float64}()
        fix_Q = Dict{Tuple{Int,Symbol},Float64}()
        for other in units
            if other.phase != ph
                fix_P[(other.node, other.phase)] = get_P_base(other)
                fix_Q[(other.node, other.phase)] = get_Q_base(other)
            end
        end

        model = build_oe_model(sys, Y, vuf_cfg; fix_P=fix_P, fix_Q=fix_Q, unit=ph_units[1])

        @variable(model, alpha >= 0)
        @objective(model, Max, alpha)

        fair_p = Dict{Int, ConstraintRef}()
        fair_q = Dict{Int, ConstraintRef}()
        for u in ph_units
            Pvar = (u.phase==:a) ? model[:P_sol_a] :
                   (u.phase==:b) ? model[:P_sol_b] : model[:P_sol_c]
            Qvar = (u.phase==:a) ? model[:Q_sol_a] :
                   (u.phase==:b) ? model[:Q_sol_b] : model[:Q_sol_c]
            fair_p[u.node] = @constraint(model,
                Pvar[u.node] - weights[u]*alpha == get_P_base(u))
            fair_q[u.node] = @constraint(model,
                Qvar[u.node]                   == get_Q_base(u))
        end

        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        for θ in angles
            dx = cos(θ)
            dy = sin(θ)

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
                @printf("    ⚠  θ=%.1f°: %s\n", rad2deg(θ), ts)
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

# ── Lexicographic max-min OE ────────────────────────────────
"""
    compute_operating_envelopes_lexmaxmin(sys, Y, vuf_cfg, baseline_res; n_angles, tol, V_base)

Pure Lexicographic Max-Min OE computation.

No sensitivity weights are used. The network physics alone determines
who binds at each stage. At each angle θ and stage s:
- All free units move by the same scalar α (no w_u scaling).
- The solver maximises α until the first network constraint binds.
- Bound units are fixed at their achieved injection; remaining
  units proceed to the next stage with a potentially larger α.

This is lex max-min optimal over the weighted step vector in the
pure sense: no external fairness parameter is imposed.
"""
function compute_operating_envelopes_lexmaxmin(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int    = 36,
    tol::Float64     = 1e-4
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units.")
    println("  Method: Pure Lexicographic Max-Min ($n_angles angles, tol=$tol)")

    get_P_base(u) = (u.phase==:a) ? baseline_res.P_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.P_sol_b[u.node] :
                                    baseline_res.P_sol_c[u.node]
    get_Q_base(u) = (u.phase==:a) ? baseline_res.Q_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.Q_sol_b[u.node] :
                                    baseline_res.Q_sol_c[u.node]

    angles  = range(0, 2π, length=n_angles+1)[1:end-1]
    oe_list = OEBoundary[]

    for ph in [:a, :b, :c]
        ph_units = filter(u -> u.phase == ph, units)
        isempty(ph_units) && continue
        n_ph = length(ph_units)

        println("  Phase $(uppercase(string(ph))) — $n_ph units")

        fix_P_base = Dict{Tuple{Int,Symbol},Float64}()
        fix_Q_base = Dict{Tuple{Int,Symbol},Float64}()
        for other in units
            if other.phase != ph
                fix_P_base[(other.node, other.phase)] = get_P_base(other)
                fix_Q_base[(other.node, other.phase)] = get_Q_base(other)
            end
        end

        ph_bounds = Dict(u => (P=Float64[], Q=Float64[]) for u in ph_units)

        for (k, θ) in enumerate(angles)
            dx = cos(θ)
            dy = sin(θ)

            allocated  = Dict{PVUnit, Float64}()
            free_units = copy(ph_units)
            warm       = Dict{Symbol, Vector{Float64}}()
            stage      = 1
            α_prev     = -Inf

            while !isempty(free_units)

                stage_fix_P = copy(fix_P_base)
                stage_fix_Q = copy(fix_Q_base)
                for (u, α_fixed) in allocated
                    stage_fix_P[(u.node, u.phase)] = get_P_base(u) + α_fixed * dx
                    stage_fix_Q[(u.node, u.phase)] = get_Q_base(u) + α_fixed * dy
                end

                model = build_oe_model(sys, Y, vuf_cfg; fix_P=stage_fix_P, fix_Q=stage_fix_Q, unit=free_units[1])

                if length(free_units) > 1
                    set_optimizer_attribute(model, "tol",      1e-4)
                    set_optimizer_attribute(model, "max_iter",  500)
                else
                    set_optimizer_attribute(model, "tol",      1e-6)
                    set_optimizer_attribute(model, "max_iter", 2000)
                end

                if !isempty(warm)
                    set_start_value.(model[:V_a_r], warm[:V_a_r])
                    set_start_value.(model[:V_a_i], warm[:V_a_i])
                    set_start_value.(model[:V_b_r], warm[:V_b_r])
                    set_start_value.(model[:V_b_i], warm[:V_b_i])
                    set_start_value.(model[:V_c_r], warm[:V_c_r])
                    set_start_value.(model[:V_c_i], warm[:V_c_i])
                    set_optimizer_attribute(model, "warm_start_init_point", "yes")
                end

                @variable(model, α_stage >= 0)
                @objective(model, Max, α_stage)

                fair_p = Dict{Int, ConstraintRef}()
                fair_q = Dict{Int, ConstraintRef}()
                for u in free_units
                    Pvar = (u.phase==:a) ? model[:P_sol_a] :
                           (u.phase==:b) ? model[:P_sol_b] : model[:P_sol_c]
                    Qvar = (u.phase==:a) ? model[:Q_sol_a] :
                           (u.phase==:b) ? model[:Q_sol_b] : model[:Q_sol_c]
                    fair_p[u.node] = @constraint(model,
                        Pvar[u.node] - α_stage * dx == get_P_base(u))
                    fair_q[u.node] = @constraint(model,
                        Qvar[u.node] - α_stage * dy == get_Q_base(u))
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

                warm[:V_a_r] = value.(model[:V_a_r])
                warm[:V_a_i] = value.(model[:V_a_i])
                warm[:V_b_r] = value.(model[:V_b_r])
                warm[:V_b_i] = value.(model[:V_b_i])
                warm[:V_c_r] = value.(model[:V_c_r])
                warm[:V_c_i] = value.(model[:V_c_i])

                α_val = value(α_stage)

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

                slacks = Dict{PVUnit, Float64}()
                for u in free_units
                    P_achieved = get_P_base(u) + α_val * dx
                    Q_achieved = get_Q_base(u) + α_val * dy
                    S_cap = (u.phase==:a) ? sys.S_Sol_Gen_ub_a[u.node] :
                            (u.phase==:b) ? sys.S_Sol_Gen_ub_b[u.node] :
                                            sys.S_Sol_Gen_ub_c[u.node]
                    slacks[u] = S_cap^2 - (P_achieved^2 + Q_achieved^2)
                end

                min_slack   = minimum(values(slacks))
                newly_bound = filter(u -> slacks[u] <= min_slack + tol, free_units)

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
                        join(["$(u.node)($(u.phase))" for u in free_units],  ","))

                stage += 1
            end

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