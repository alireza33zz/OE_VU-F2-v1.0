# ============================================================
# 3-PHASE AC-OPF MODEL (rectangular voltages, exact power balance)
# ============================================================

"""
    add_vuf_constraints!(model, vuf_cfg, V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i)

Inline VUF constraints without auxiliary variables:
  enforce |V₋|² ≤ ε²·|V₊|² and (|V₊|² ≥ floor) per bus.
"""
function add_vuf_constraints!(
    model::JuMP.Model,
    vuf_cfg::VUFConfig,
    V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i
)
    isempty(vuf_cfg.buses) && return
    ar, ai   = -0.5,  sqrt(3)/2
    a2r, a2i = -0.5, -sqrt(3)/2
    ε2  = vuf_cfg.threshold^2
    flr = vuf_cfg.pos_seq_floor

    for b in vuf_cfg.buses
        # Positive-seq
        Vp_r = @expression(model, (1/3)*( V_a_r[b] + ar*V_b_r[b] - ai*V_b_i[b] + a2r*V_c_r[b] - a2i*V_c_i[b] ))
        Vp_i = @expression(model, (1/3)*( V_a_i[b] + ar*V_b_i[b] + ai*V_b_r[b] + a2r*V_c_i[b] + a2i*V_c_r[b] ))
        # Negative-seq
        Vn_r = @expression(model, (1/3)*( V_a_r[b] + a2r*V_b_r[b] - a2i*V_b_i[b] + ar*V_c_r[b] - ai*V_c_i[b] ))
        Vn_i = @expression(model, (1/3)*( V_a_i[b] + a2r*V_b_i[b] + a2i*V_b_r[b] + ar*V_c_i[b] + ai*V_c_r[b] ))

        # Floor on |V+|^2 and the VUF inequality
        @constraint(model, Vp_r^2 + Vp_i^2 >= flr)
        @constraint(model, Vn_r^2 + Vn_i^2 <= ε2 * (Vp_r^2 + Vp_i^2))
    end
end

"""
    add_power_balance_exact!(model, sys, Y, ...)

Exact complex power balance  S = diag(V) * conj(Y * V)
with V in rectangular coordinates, built with JuMP-safe expressions.
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

    # Split Y into real G and imag B
    G = real.(Y)
    B = imag.(Y)

    # Helpers to fetch the VariableRef by stacked index
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

    # J = Y * conj(V) = (G + jB) * (V_r - j V_i)
    #   J_r[k] = Σ_j ( G[k,j]*V_r[j] + B[k,j]*V_i[j] )
    #   J_i[k] = Σ_j ( B[k,j]*V_r[j] - G[k,j]*V_i[j] )
    @expression(model, J_r[k=1:n3], sum(G[k,j]*vr(j) + B[k,j]*vi(j) for j in 1:n3))
    @expression(model, J_i[k=1:n3], sum(B[k,j]*vr(j) - G[k,j]*vi(j) for j in 1:n3))

    # S = V ⊙ J  (since conj was applied inside J):
    #   P[k] =  V_r[k]*J_r[k] - V_i[k]*J_i[k]
    #   Q[k] =  V_r[k]*J_i[k] + V_i[k]*J_r[k]
    @expression(model, P[k=1:n3], vr(k)*J_r[k] - vi(k)*J_i[k])
    @expression(model, Q[k=1:n3], vr(k)*J_i[k] + vi(k)*J_r[k])

    # Map back to (node, phase) and impose balances
    for i in 1:nN
        ka, kb, kc = 3i-2, 3i-1, 3i

        @constraint(model, P_G_a[i] + P_sol_a[i] - sys.P_Demand_a[i] == P[ka])
        @constraint(model, P_G_b[i] + P_sol_b[i] - sys.P_Demand_b[i] == P[kb])
        @constraint(model, P_G_c[i] + P_sol_c[i] - sys.P_Demand_c[i] == P[kc])

        @constraint(model, Q_G_a[i] + Q_sol_a[i] - sys.Q_Demand_a[i] == Q[ka])
        @constraint(model, Q_G_b[i] + Q_sol_b[i] - sys.Q_Demand_b[i] == Q[kb])
        @constraint(model, Q_G_c[i] + Q_sol_c[i] - sys.Q_Demand_c[i] == Q[kc])
    end
end

"""
    build_base_model(sys, Y, vuf_cfg; fix_P, fix_Q, max_iter, tol) -> JuMP.Model

Variables + all network constraints, no objective, not solved.
Shared by the baseline OPF and every OE sub-problem.

`fix_P` / `fix_Q` map (node, phase) → fixed PV injection.
"""
function build_base_model(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig = VUFConfig();
    fix_P::Dict{Tuple{Int,Symbol},Float64} = Dict{Tuple{Int,Symbol},Float64}(),
    fix_Q::Dict{Tuple{Int,Symbol},Float64} = Dict{Tuple{Int,Symbol},Float64}(),
    max_iter::Int = 200_000,
    tol::Float64  = 1e-6
)::JuMP.Model

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", max_iter)
    set_optimizer_attribute(model, "tol",      tol)
    set_optimizer_attribute(model, "print_level", 0)
    set_optimizer_attribute(model, "linear_solver", "mumps")
    set_optimizer_attribute(model, "warm_start_init_point", "no")

    nN = sys.nN
    V0 = sys.V_nom

    # ----------------------------------------------------------
    # Variables
    # ----------------------------------------------------------
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

    # ----------------------------------------------------------
    # Slack bus — fixed balanced reference at node 1
    # ----------------------------------------------------------
    fix(V_a_r[1], V0*cos(0);     force=true)
    fix(V_a_i[1], V0*sin(0);     force=true)
    fix(V_b_r[1], V0*cos(-2π/3); force=true)
    fix(V_b_i[1], V0*sin(-2π/3); force=true)
    fix(V_c_r[1], V0*cos( 2π/3); force=true)
    fix(V_c_i[1], V0*sin( 2π/3); force=true)

    # ----------------------------------------------------------
    # Conventional generator limits
    # ----------------------------------------------------------
    for i in 1:nN
        s = sys.Gen_Status[i]
        for (PG, QG) in ((P_G_a,Q_G_a),(P_G_b,Q_G_b),(P_G_c,Q_G_c))
            @constraint(model, sys.P_Gen_lb[i]*s <= PG[i] <= sys.P_Gen_ub[i]*s)
            @constraint(model, sys.Q_Gen_lb[i]*s <= QG[i] <= sys.Q_Gen_ub[i]*s)
        end
    end

    # ----------------------------------------------------------
    # PV apparent power limits  |S_sol|² ≤ S_max²
    # ----------------------------------------------------------
    @constraint(model,[i in 1:nN], P_sol_a[i]^2+Q_sol_a[i]^2 <= sys.S_Sol_Gen_ub_a[i]^2)
    @constraint(model,[i in 1:nN], P_sol_b[i]^2+Q_sol_b[i]^2 <= sys.S_Sol_Gen_ub_b[i]^2)
    @constraint(model,[i in 1:nN], P_sol_c[i]^2+Q_sol_c[i]^2 <= sys.S_Sol_Gen_ub_c[i]^2)

    # ----------------------------------------------------------
    # Fixed PV setpoints (OE scans: freeze units that are not being traced)
    # ----------------------------------------------------------
    for ((node, ph), val) in fix_P
        fix(_pick(ph, P_sol_a, P_sol_b, P_sol_c)[node], val; force=true)
    end
    for ((node, ph), val) in fix_Q
        fix(_pick(ph, Q_sol_a, Q_sol_b, Q_sol_c)[node], val; force=true)
    end

    # ----------------------------------------------------------
    # Voltage magnitude limits  Vmin² ≤ |V|² ≤ Vmax²
    # ----------------------------------------------------------
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_a_r[i]^2+V_a_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_b_r[i]^2+V_b_i[i]^2 <= sys.V_ub[i]^2)
    @constraint(model,[i in 1:nN], sys.V_lb[i]^2 <= V_c_r[i]^2+V_c_i[i]^2 <= sys.V_ub[i]^2)

    # ----------------------------------------------------------
    # VUF constraints (optional, per-bus) and exact nodal power balance
    # ----------------------------------------------------------
    add_vuf_constraints!(model, vuf_cfg, V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i)
    add_power_balance_exact!(model, sys, Y,
        V_a_r, V_a_i, V_b_r, V_b_i, V_c_r, V_c_i,
        P_G_a, P_G_b, P_G_c, Q_G_a, Q_G_b, Q_G_c,
        P_sol_a, P_sol_b, P_sol_c, Q_sol_a, Q_sol_b, Q_sol_c)

    return model
end

"""
    build_and_solve_opf(sys, Y, vuf_cfg; max_iter, tol) -> JuMP.Model

Baseline OPF: minimise conventional generation cost.
"""
function build_and_solve_opf(
    sys::SystemData,
    Y::Matrix{ComplexF64},
    vuf_cfg::VUFConfig = VUFConfig();
    max_iter::Int = 200_000,
    tol::Float64  = 1e-6
)::JuMP.Model

    model = build_base_model(sys, Y, vuf_cfg; max_iter=max_iter, tol=tol)
    P_G_a, P_G_b, P_G_c = model[:P_G_a], model[:P_G_b], model[:P_G_c]
    @objective(model, Min,
        sum(sys.P_Cost_G[i]*(P_G_a[i]+P_G_b[i]+P_G_c[i]) for i in 1:sys.nN))

    optimize!(model)
    return model
end

"""
    extract_result(model, nN) -> OPFResult
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

# ============================================================
# POST-PROCESSING: VUF & LINE FLOWS
# ============================================================

"""
    compute_vuf(res, nN) -> Vector{Float64}

Voltage Unbalance Factor (%) at every bus:  VUF(%) = 100 · |V₋| / |V₊|
"""
function compute_vuf(res::OPFResult, nN::Int)::Vector{Float64}
    ar, ai   = -0.5,  sqrt(3)/2
    a2r, a2i = -0.5, -sqrt(3)/2

    vuf = zeros(nN)
    for i in 1:nN
        Va_r = res.V_r_a[i];  Va_i = res.V_i_a[i]
        Vb_r = res.V_r_b[i];  Vb_i = res.V_i_b[i]
        Vc_r = res.V_r_c[i];  Vc_i = res.V_i_c[i]

        # Positive sequence
        Vp_r = (1/3)*(Va_r + ar*Vb_r - ai*Vb_i + a2r*Vc_r - a2i*Vc_i)
        Vp_i = (1/3)*(Va_i + ar*Vb_i + ai*Vb_r + a2r*Vc_i + a2i*Vc_r)

        # Negative sequence
        Vn_r = (1/3)*(Va_r + a2r*Vb_r - a2i*Vb_i + ar*Vc_r - ai*Vc_i)
        Vn_i = (1/3)*(Va_i + a2r*Vb_i + a2i*Vb_r + ar*Vc_i + ai*Vc_r)

        Vp_mag = sqrt(Vp_r^2 + Vp_i^2)
        Vn_mag = sqrt(Vn_r^2 + Vn_i^2)

        vuf[i] = Vp_mag > 1e-9 ? 100.0 * Vn_mag / Vp_mag : 0.0
    end
    return vuf
end

"""
    compute_line_flows(sys, res) -> (Pij, Pji, Qij, Qji)

Per-phase sending/receiving-end flows (nL × 3) using the π-model.
"""
function compute_line_flows(sys::SystemData, res::OPFResult)
    nL = sys.nL

    Pij = zeros(nL, 3);  Pji = zeros(nL, 3)
    Qij = zeros(nL, 3);  Qji = zeros(nL, 3)

    for k in 1:nL
        # 6×6 branch admittance (series admittance + half shunt each end)
        Yl  = inv(sys.Zij[k])
        Yl6 = [(Yl + 0.5 .* sys.Ysh) (-Yl);
               (-Yl)                (Yl + 0.5 .* sys.Ysh)]

        i = sys.fbus[k]
        j = sys.tbus[k]

        # [Vi_abc; Vj_abc]
        v = [res.V_r_a[i] + im*res.V_i_a[i];
             res.V_r_b[i] + im*res.V_i_b[i];
             res.V_r_c[i] + im*res.V_i_c[i];
             res.V_r_a[j] + im*res.V_i_a[j];
             res.V_r_b[j] + im*res.V_i_b[j];
             res.V_r_c[j] + im*res.V_i_c[j]]

        # S = V * conj(I), where I = Y * V
        Sij = v[1:3] .* conj.(Yl6[1:3, :] * v)
        Sji = v[4:6] .* conj.(Yl6[4:6, :] * v)

        Pij[k, :] = real.(Sij);  Qij[k, :] = imag.(Sij)
        Pji[k, :] = real.(Sji);  Qji[k, :] = imag.(Sji)
    end

    return Pij, Pji, Qij, Qji
end
