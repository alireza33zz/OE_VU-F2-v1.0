# ============================================================
# RESULT REPORTING
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
