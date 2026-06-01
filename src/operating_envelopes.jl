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

    # Build the full baseline voltage vector (3N) in node-interleaved order
    # index 3*(i-1)+1 = phase a of node i
    # index 3*(i-1)+2 = phase b of node i
    # index 3*(i-1)+3 = phase c of node i
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

    # Phase-local normalisation — unchanged
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

# ============================================================
# 10. OE PLOTS
# ============================================================

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

        # Hardware S_max circle
        plot!(ax, S_max.*cos.(θ_circ), S_max.*sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label="S_max")

        # OE shaded polygon
        if length(oe.P_boundary) >= 3
            cx = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
            idx = sortperm(atan.(oe.Q_boundary.-cy, oe.P_boundary.-cx))
            Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
            Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale
            plot!(ax, Ps, Qs; seriestype=:shape,
                  fillcolor=:dodgerblue, fillalpha=0.18,
                  lc=:dodgerblue, lw=1.8, label="OE")
        end

        # Baseline operating point
        scatter!(ax, [oe.P_baseline*scale], [oe.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label="Baseline")

        # Zero-crossing lines
        hline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)
        vline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)

        title!(ax, "Node $(oe.unit.node) · Ph $(uppercase(string(oe.unit.phase)))")
        xlabel!(ax, "P (kW)");  ylabel!(ax, "Q (kvar)")
        plot!(ax, legend=(k==1), legendfont=font(7), titlefont=font(9,:bold),
              tickfont=font(7), guidefont=font(8),
              aspect_ratio=:equal, framestyle=:box, grid=true, gridalpha=0.25)
    end

# ── Weighted-fairness OE ────────────────────────────────────
function compute_operating_envelopes(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig,
    baseline_res::OPFResult;
    n_angles::Int    = 36,
    V_base::Float64  = 400.0 / sqrt(3)
)::Vector{OEBoundary}

    units   = find_pv_units(sys)
    n_units = length(units)
    println("  Found $n_units PV units. Tracing OE ($n_angles angles, sensitivity-weighted fairness)...")

    get_P_base(u) = (u.phase==:a) ? baseline_res.P_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.P_sol_b[u.node] : baseline_res.P_sol_c[u.node]
    get_Q_base(u) = (u.phase==:a) ? baseline_res.Q_sol_a[u.node] :
                    (u.phase==:b) ? baseline_res.Q_sol_b[u.node] : baseline_res.Q_sol_c[u.node]

    # Pre-compute Z = Y^{-1} once; pass cached matrix to weight function
    Z_full = build_impedance_matrix(Y; slack_node=1)

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

        # Fix every unit on other phases at their baseline
        fix_P = Dict{Tuple{Int,Symbol},Float64}()
        fix_Q = Dict{Tuple{Int,Symbol},Float64}()
        for other in units
            if other.phase != ph
                fix_P[(other.node, other.phase)] = get_P_base(other)
                fix_Q[(other.node, other.phase)] = get_Q_base(other)
            end
        end

        model = build_oe_model(sys, Y, vuf_cfg;
                               fix_P=fix_P, fix_Q=fix_Q, unit=ph_units[1])

        @variable(model, alpha >= 0)
        @objective(model, Max, alpha)

        # Fairness constraints: unit i steps by w_i * alpha in the shared direction (cos θ, sin θ).
        # Coefficients on alpha are updated each angle via set_normalized_coefficient.
        fair_p = Dict{Int, ConstraintRef}()
        fair_q = Dict{Int, ConstraintRef}()
        for u in ph_units
            Pvar = (u.phase==:a) ? model[:P_sol_a] :
                   (u.phase==:b) ? model[:P_sol_b] : model[:P_sol_c]
            Qvar = (u.phase==:a) ? model[:Q_sol_a] :
                   (u.phase==:b) ? model[:Q_sol_b] : model[:Q_sol_c]
            # Initial coefficients (θ=0): P = P_base + w*alpha, Q = Q_base
            fair_p[u.node] = @constraint(model,
                Pvar[u.node] - weights[u]*alpha == get_P_base(u))
            fair_q[u.node] = @constraint(model,
                Qvar[u.node]                   == get_Q_base(u))
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

"""
    compute_operating_envelopes_lexmaxmin(sys, Y, vuf_cfg, baseline_res; n_angles, tol)

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
    n_angles::Int  = 36,
    tol::Float64   = 1e-4
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

        # Fix all other-phase units at baseline — unchanged from before
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

            # allocated[(u)] = the α value at which unit u was fixed
            # In pure lex: the step of unit u is exactly α (no w_u)
            allocated = Dict{PVUnit, Float64}()
            free_units = copy(ph_units)

            # ADD: warm start storage, initialised empty for first stage
            warm = Dict{Symbol, Vector{Float64}}()

            stage = 1
            α_prev = -Inf
            while !isempty(free_units)

                # Build fix dicts: baseline fixes for other phases
                # + achieved injection fixes for already-bound units
                stage_fix_P = copy(fix_P_base)
                stage_fix_Q = copy(fix_Q_base)
                for (u, α_fixed) in allocated
                    # Bound unit injection is fixed at its achieved value
                    stage_fix_P[(u.node, u.phase)] =
                        get_P_base(u) + α_fixed * dx
                    stage_fix_Q[(u.node, u.phase)] =
                        get_Q_base(u) + α_fixed * dy
                end

                model = build_oe_model(sys, Y, vuf_cfg;
                                       fix_P = stage_fix_P,
                                       fix_Q = stage_fix_Q,
                                       unit  = free_units[1])

# ADD: relax tolerance for intermediate stages
if length(free_units) > 1
    set_optimizer_attribute(model, "tol",      1e-4)
    set_optimizer_attribute(model, "max_iter",  500)
else
    # Final stage: full precision for accurate boundary point
    set_optimizer_attribute(model, "tol",      1e-6)
    set_optimizer_attribute(model, "max_iter", 2000)
end

        # ADD: apply warm start if available from previous stage
        if !isempty(warm)
            set_start_value.(model[:V_a_r], warm[:V_a_r])
            set_start_value.(model[:V_a_i], warm[:V_a_i])
            set_start_value.(model[:V_b_r], warm[:V_b_r])
            set_start_value.(model[:V_b_i], warm[:V_b_i])
            set_start_value.(model[:V_c_r], warm[:V_c_r])
            set_start_value.(model[:V_c_i], warm[:V_c_i])
            set_optimizer_attribute(model, "warm_start_init_point", "yes")
        end                       

                # Shared stage scalar — NO weight scaling
                @variable(model, α_stage >= 0)
                @objective(model, Max, α_stage)

                # Every free unit moves by exactly α in direction (dx, dy)
                fair_p = Dict{Int, ConstraintRef}()
                fair_q = Dict{Int, ConstraintRef}()
                for u in free_units
                    Pvar = (u.phase==:a) ? model[:P_sol_a] :
                           (u.phase==:b) ? model[:P_sol_b] : model[:P_sol_c]
                    Qvar = (u.phase==:a) ? model[:Q_sol_a] :
                           (u.phase==:b) ? model[:Q_sol_b] : model[:Q_sol_c]

                    # P_u = P_base_u + α * dx
                    fair_p[u.node] = @constraint(model,
                        Pvar[u.node] - α_stage * dx == get_P_base(u))
                    # Q_u = Q_base_u + α * dy
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

                # ADD: store solution as warm start for next stage
                warm[:V_a_r] = value.(model[:V_a_r])
                warm[:V_a_i] = value.(model[:V_a_i])
                warm[:V_b_r] = value.(model[:V_b_r])
                warm[:V_b_i] = value.(model[:V_b_i])
                warm[:V_c_r] = value.(model[:V_c_r])
                warm[:V_c_i] = value.(model[:V_c_i])

                α_val = value(α_stage)

                # ADD THIS BLOCK immediately after α_val = value(α_stage):
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

                # Determine which free units are binding at this stage.
                # In pure lex all free units share the same α, so the
                # binding condition is purely whether the network constraint
                # active at α_val originates from this unit.
                # We detect this by checking which units would violate a
                # network constraint if α were increased infinitesimally —
                # in practice we simply fix ALL free units at this stage's
                # α_val and move to the next stage, because the solver has
                # already found the maximum feasible α for the entire
                # free set simultaneously.
                #
                # The lex structure then comes from the fact that in the
                # next stage the bound units are hard-fixed, so the
                # remaining free units can potentially achieve a larger α.

                # Identify units that are at their individual capacity limit:
                # query each unit's solar/voltage slack to find who is tight.
                # Simpler and equally correct: bind the unit(s) with the
                # smallest individual headroom, measured by checking whether
                # fixing them and re-solving gives a strictly larger α.
                #
                # For robustness we use the tolerance-based approach:
                # all free units receive α_val at this stage; any unit
                # whose capacity constraint is active (detected via dual or
                # by a single re-solve probe) is bound. In the absence of
                # dual access through Ipopt, we bind the unit with the
                # minimum slack in its apparent-power constraint.

                # Compute apparent-power slack for each free unit
                slacks = Dict{PVUnit, Float64}()
                for u in free_units
                    P_achieved = get_P_base(u) + α_val * dx
                    Q_achieved = get_Q_base(u) + α_val * dy
                    S_cap = (u.phase==:a) ? sys.S_Sol_Gen_ub_a[u.node] :
                            (u.phase==:b) ? sys.S_Sol_Gen_ub_b[u.node] :
                                            sys.S_Sol_Gen_ub_c[u.node]
                    slacks[u] = S_cap^2 - (P_achieved^2 + Q_achieved^2)
                end

                min_slack  = minimum(values(slacks))
                newly_bound = filter(u -> slacks[u] <= min_slack + tol, free_units)

                # Guard: always bind at least one unit per stage
                if isempty(newly_bound)
                    newly_bound = [free_units[argmin([slacks[u]
                                   for u in free_units])]]
                end

                for u in newly_bound
                    allocated[u] = α_val
                end
                filter!(u -> u ∉ newly_bound, free_units)

                @printf("    θ=%5.1f° stage %d: α*=%.5f | bound: %s | free: %s\n",
                        rad2deg(θ), stage, α_val,
                        join(["$(u.node)($(u.phase))"
                              for u in newly_bound], ","),
                        join(["$(u.node)($(u.phase))"
                              for u in free_units], ","))

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

# Signed shortest angular difference  b - a  in (-π, π]
angular_diff(b, a) = mod(b - a + π, 2π) - π


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

    # Build the full baseline voltage vector (3N) in node-interleaved order
    # index 3*(i-1)+1 = phase a of node i
    # index 3*(i-1)+2 = phase b of node i
    # index 3*(i-1)+3 = phase c of node i
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

    # Phase-local normalisation — unchanged
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

# ============================================================
# 10. OE PLOTS
# ============================================================

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

        # Hardware S_max circle
        plot!(ax, S_max.*cos.(θ_circ), S_max.*sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label="S_max")

        # OE shaded polygon
        if length(oe.P_boundary) >= 3
            cx = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
            idx = sortperm(atan.(oe.Q_boundary.-cy, oe.P_boundary.-cx))
            Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
            Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale
            plot!(ax, Ps, Qs; seriestype=:shape,
                  fillcolor=:dodgerblue, fillalpha=0.18,
                  lc=:dodgerblue, lw=1.8, label="OE")
        end

        # Baseline operating point
        scatter!(ax, [oe.P_baseline*scale], [oe.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label="Baseline")

        # Zero-crossing lines
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

# ============================================================
# 11. RESULT REPORTING
# ============================================================

"""
    print_results(sys, res, vuf_buses)

Print all result tables to stdout.
"""
function print_results(sys::SystemData, res::OPFResult, vuf_buses::Vector{Int}=Int[])
    println("\n", "="^72)
    println("  3-PHASE AC-OPF RESULTS")
    println("="^72)
    println("  Status    : ", res.status)
    println("  Objective : ", round(res.obj; digits=2), "  (total generation cost)")
    println("="^72, "\n")

    nN  = sys.nN

    # ---- Voltage magnitudes & angles ----
    mag(r,i) = sqrt.(r.^2 .+ i.^2)
    arg(r,i) = rad2deg.(atan.(i, r))

    Va_mag = mag(res.V_r_a, res.V_i_a);  Va_arg = arg(res.V_r_a, res.V_i_a)
    Vb_mag = mag(res.V_r_b, res.V_i_b);  Vb_arg = arg(res.V_r_b, res.V_i_b)
    Vc_mag = mag(res.V_r_c, res.V_i_c);  Vc_arg = arg(res.V_r_c, res.V_i_c)

    # ---- VUF (post-processing) ----
    vuf = compute_vuf(res, nN)

    # Flag buses where VUF is constrained
    flag = [i in vuf_buses ? "*" : " " for i in 1:nN]

    println("── Nodal voltages — magnitude (V), angle (°), VUF (%)  [* = VUF constrained] ──")
    println(DataFrame(
        Node    = 1:nN,
        Flag    = flag,
        Va_mag  = round.(Va_mag; digits=2),
        Va_arg  = round.(Va_arg; digits=1),
        Vb_mag  = round.(Vb_mag; digits=2),
        Vb_arg  = round.(Vb_arg; digits=1),
        Vc_mag  = round.(Vc_mag; digits=2),
        Vc_arg  = round.(Vc_arg; digits=1),
        VUF_pct = round.(vuf;    digits=2)
    ))

    println("\n── Conventional generator dispatch (W / var) ──")
    println(DataFrame(
        Node = 1:nN,
        PG_a = round.(res.P_G_a; digits=2),  QG_a = round.(res.Q_G_a; digits=2),
        PG_b = round.(res.P_G_b; digits=2),  QG_b = round.(res.Q_G_b; digits=2),
        PG_c = round.(res.P_G_c; digits=2),  QG_c = round.(res.Q_G_c; digits=2)))

    println("\n── PV generator dispatch (W / var) ──")
    println(DataFrame(
        Node   = 1:nN,
        Psol_a = round.(res.P_sol_a; digits=2),  Qsol_a = round.(res.Q_sol_a; digits=2),
        Psol_b = round.(res.P_sol_b; digits=2),  Qsol_b = round.(res.Q_sol_b; digits=2),
        Psol_c = round.(res.P_sol_c; digits=2),  Qsol_c = round.(res.Q_sol_c; digits=2)))

    Pij, Pji, Qij, Qji = compute_line_flows(sys, res)

    println("\n── Active power flows (W)  [ij = sending end, ji = receiving end] ──")
    println(DataFrame(
        Line  = 1:sys.nL, From = sys.fbus, To = sys.tbus,
        Pij_a = round.(Pij[:,1]; digits=2), Pij_b = round.(Pij[:,2]; digits=2),
        Pij_c = round.(Pij[:,3]; digits=2),
        Pji_a = round.(Pji[:,1]; digits=2), Pji_b = round.(Pji[:,2]; digits=2),
        Pji_c = round.(Pji[:,3]; digits=2)))

    println("\n── Reactive power flows (var) ──")
    println(DataFrame(
        Line  = 1:sys.nL, From = sys.fbus, To = sys.tbus,
        Qij_a = round.(Qij[:,1]; digits=2), Qij_b = round.(Qij[:,2]; digits=2),
        Qij_c = round.(Qij[:,3]; digits=2),
        Qji_a = round.(Qji[:,1]; digits=2), Qji_b = round.(Qji[:,2]; digits=2),
        Qji_c = round.(Qji[:,3]; digits=2)))
end

# ============================================================
# ADD-ON: BASELINE OPF DASHBOARD  (4 figures)  —  WHITE THEME
# ============================================================
# Paste at the end of opf_3phase_with_OE_VUF.jl
# Requires: Plots.jl + Printf (both already available)
#
# Call after solving the baseline:
#   figs = plot_dashboard(sys, baseline_res)
#   for (name, f) in figs; display(f); end
# ============================================================

using Printf

# ─────────────────────────────────────────────────────────────────────────────
# SHARED STYLE CONSTANTS  —  WHITE THEME
# ─────────────────────────────────────────────────────────────────────────────

# White canvas
const _BG       = RGB(1, 1, 1)               # figure canvas
const _BG2      = RGB(1, 1, 1)               # panel surface
const _GRID     = RGBA(0,0,0,0.12)           # light grey grid
const _FG       = RGB(0.10, 0.12, 0.16)      # primary text (dark)
const _FG2      = RGB(0.35, 0.37, 0.43)      # secondary / axis text

# Phase colours (saturated, work on white bg)
const _PA  = RGB(0.937, 0.380, 0.216)        # #EF6137  coral-red   phase A
const _PB  = RGB(0.149, 0.737, 0.557)        # #26BC8E  mint-teal   phase B
const _PC  = RGB(0.290, 0.604, 0.976)        # #4A9AF9  sky-blue    phase C

# Pale versions (used for demand/PV fills)
const _PA2 = RGBA(0.937, 0.380, 0.216, 0.28)
const _PB2 = RGBA(0.149, 0.737, 0.557, 0.28)
const _PC2 = RGBA(0.290, 0.604, 0.976, 0.28)

# Demand (strong, saturated)
const _PA_D = RGB(0.87, 0.30, 0.09)   # deep orange (Phase A)
const _PB_D = RGB(0.05, 0.60, 0.45)   # teal-green  (Phase B)
const _PC_D = RGB(0.16, 0.47, 0.86)   # sky-blue    (Phase C)

# PV fills (near-white tints of the same hues)
const _PA_PV_FILL = RGB(0.99, 0.94, 0.90)  # very light orange
const _PB_PV_FILL = RGB(0.90, 0.97, 0.94)  # very light teal
const _PC_PV_FILL = RGB(0.92, 0.96, 1.00)  # very light blue

# Optional: common hatch for PV (works best on GR backend)
const _PV_HATCH = :x

# Voltage limit band colour (subtle on white)
const _VLIM = RGBA(0.98, 0.58, 0.18, 0.2)

# VUF accent
const _VUF_HI = RGB(0.98, 0.05, 0.05)        # orange — VUF > threshold

# Common figure kwargs
_FIG_KW(W, H) = (
    size                    = (W, H),
    background_color        = _BG,
    background_color_inside = _BG2,
    foreground_color        = _FG,
    foreground_color_axis   = _FG2,
    foreground_color_border = _FG2,
    foreground_color_grid   = _GRID,
    foreground_color_text   = _FG,
    foreground_color_guide  = _FG,
    dpi                     = 160,
    grid                    = true,
    gridalpha               = 0.14,
    gridlinewidth           = 0.6,
    minorgrid               = false,
    framestyle              = :box,
    tick_direction          = :out,
    fontfamily              = "DejaVu Sans",
)

# Thin spine helper: draw a horizontal zero-line (dark on white)
_zeroline!(p) = hline!(p, [0.0]; lc=RGBA(0,0,0,0.25), lw=0.8, label=false)

# Node x-positions  (1..18, integer ticks)
const _NODES = 1:13


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — ACTIVE POWER (P demand ↓  /  P injection ↑)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_active_power(sys, res) -> Plots.Plot

Three-column cluster per node (Ph-A, Ph-B, Ph-C). Within each phase, two bars:
  ▼  Active demand  (negative → downward)
  ▲  PV injection   (positive → upward)

x-axis: node index.  y-axis: kW.
Legend is arranged as 3 columns × 2 rows and placed outside (top center).
Solar injection uses a distinct hatch fill.
"""
function _fig_active_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3   # W → kW

    # --- build bar matrices (nN × 3) -----------------------------------------
    Pd = hcat(sys.P_Demand_a, sys.P_Demand_b, sys.P_Demand_c) .* scale
    Pp = hcat(res.P_sol_a,    res.P_sol_b,    res.P_sol_c)    .* scale

    # Clamp near-zero to exactly 0 for clean rendering
    Pd[abs.(Pd) .< 1e-2] .= 0.0
    Pp[abs.(Pp) .< 1e-2] .= 0.0

    cols_solid = [_PA, _PB, _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Demand Ph-A", "Demand Ph-B", "Demand Ph-C"]
    labels_pv  = ["PV Ph-A",     "PV Ph-B",     "PV Ph-C"]

    # Phase cluster offsets (A,B,C) → (−δ, 0, +δ)
    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)

    # Inner separation within a phase (Demand vs PV)
    δ_inner   = 0.09
    w_bar     = 0.16

    p = plot(; _FIG_KW(1320, 520)...,
        title        = "Active Power — Demand (↓) and PV Injection (↑)",
        titlefont    = (11, _FG, :bold),
        xlabel       = "Node",
        ylabel       = "P  [kW]",
        xticks       = (nodes, string.(nodes)),
        xrotation    = 0,
        ylims        = (-maximum(abs.(Pd))*1.1, maximum(Pp)*1.1),
        # Legend requested outside, top center
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legendfontcolor = _FG2,
        legend_columns = 6,   # 3×2 layout for 6 entries
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    # Demand/PV bars grouped by phase with small inner offset
    for ph in 1:3
        # Pick colors per phase
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        # Demand (solid, saturated)
        bar!(p, x_base .- δ_inner, -Pd[:, ph];   # or -Qd for Fig. 2
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            fillalpha = 1.0,
            label     = labels_dem[ph])

        # PV (light tint + hatch + colored outline)
        bar!(p, x_base .+ δ_inner,  Pp[:, ph];   # or  Qp for Fig. 2
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,    # outline matches the phase’s demand color
            linewidth   = 1.1,
            fillalpha   = 1.0,
            fillpattern = _PV_HATCH,     # makes PV immediately distinct
            label       = labels_pv[ph])
    end

    return p
end


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — REACTIVE POWER (Q demand ↓  /  Q injection ↑)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_reactive_power(sys, res) -> Plots.Plot

Same layout as Fig 1, but for reactive power (kvar).
Three-column cluster per node; legend outside (top center), 3×2 layout.
Solar injection uses a distinct hatch fill.
"""
function _fig_reactive_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3

    Qd = hcat(sys.Q_Demand_a, sys.Q_Demand_b, sys.Q_Demand_c) .* scale
    Qp = hcat(res.Q_sol_a,    res.Q_sol_b,    res.Q_sol_c)    .* scale

    Qd[abs.(Qd) .< 1e-4] .= 0.0
    Qp[abs.(Qp) .< 1e-4] .= 0.0

    cols_solid = [_PA, _PB, _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Q-demand Ph-A", "Q-demand Ph-B", "Q-demand Ph-C"]
    labels_pv  = ["Q-PV Ph-A",     "Q-PV Ph-B",     "Q-PV Ph-C"]

    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    δ_inner   = 0.09
    w_bar     = 0.16

    p = plot(; _FIG_KW(1320, 520)...,
        title        = "Reactive Power — Demand (↓) and PV Injection (↑)",
        titlefont    = (11, _FG, :bold),
        xlabel       = "Node",
        ylabel       = "Q  [kvar]",
        xticks       = (nodes, string.(nodes)),
        ylims        = (-maximum(abs.(Qd))*1.1, maximum(Qp)*1.1),
        # Legend requested outside, top center
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legendfontcolor = _FG2,
        legend_columns = 6,
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    for ph in 1:3
        # Phase-based colors (same spectrum; demand strong, PV very light)
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        # Demand (downward, saturated fill)
        bar!(p, x_base .- δ_inner, -Qd[:, ph];
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            fillalpha = 1.0,
            label     = labels_dem[ph])

        # PV (upward, very light tint + colored outline + hatch)
        bar!(p, x_base .+ δ_inner,  Qp[:, ph];
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,   # phase-matched outline
            linewidth   = 1.1,
            fillalpha   = 1.0,
            fillpattern = _PV_HATCH,    # e.g., :x (requires GR or a supported backend)
            label       = labels_pv[ph])
    end

    return p
end

# ── Lexicographic max-min OE ────────────────────────────────
"""
    compute_operating_envelopes_lexmaxmin(sys, Y, vuf_cfg, baseline_res; n_angles, tol)

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
    n_angles::Int  = 36,
    tol::Float64   = 1e-4
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
