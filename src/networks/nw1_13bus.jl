# ============================================================
# NW1 — 13-bus LV feeder (400 V line-to-line)
# ============================================================

function build_system_data_13bus()::SystemData
    nN = 13;  nL = 12

    # Demand_scaling kept at 1.0 — unbalance comes from phase asymmetry, not raw
    # magnitude. Higher scaling risks voltage violations before VUF reaches the target.
    Demand_scaling = 1.0

    # Heavily loading Phase A, moderating Phase B, and stripping Phase C.
    # This creates a significant neutral point shift to hit that 1-2% VUF target.
    
    P_Demand_a = Float64.([ 0,    0,    0,    0, 3000,    0,    0, 3000,    0,    0, 3000,    0,    0]) * Demand_scaling
    P_Demand_b = Float64.([ 0,    0,    0,    0,    0, 2000,    0,    0, 2000,    0,    0, 2000,    0]) * Demand_scaling
    P_Demand_c = Float64.([ 0,    0,    0,    0,    0,    0, 1000,    0,    0, 1000,    0,    0, 1000]) * Demand_scaling
    # Maintaining Q ≈ 0.25 * P to keep power factor consistent (~0.97)
    Q_Demand_a = P_Demand_a * 0.25
    Q_Demand_b = P_Demand_b * 0.25
    Q_Demand_c = P_Demand_c * 0.25

    # Solar generation remains untouched as requested
    S_Sol_Gen_ub_a = Float64.([ 0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500,    0,    0]) .* 1.0
    S_Sol_Gen_ub_b = Float64.([ 0,    0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500,    0]) .* 1.0
    S_Sol_Gen_ub_c = Float64.([ 0,    0,    0,    0,    0,    0, 2500,    0,    0, 2500,    0,    0, 2500]) .* 1.0

    P_Cost_G   = [1.0; zeros(nN-1)]
    P_Gen_lb   = [-3_000_000.0; zeros(nN-1)]
    P_Gen_ub   = [ 3_000_000.0; zeros(nN-1)]
    Q_Gen_lb   = [-500_000.0; zeros(nN-1)]
    Q_Gen_ub   = [ 500_000.0; zeros(nN-1)]
    Gen_Status = [1; zeros(Int,nN-1)]

    V_base = 400.0 / sqrt(3)
    V_lb   = fill(0.94 * V_base, nN) # 217.1 V
    V_ub   = fill(1.10 * V_base, nN) # 254.0 V

    fbus = [1, 2, 2, 2, 2, 3, 3, 3, 3, 4,  4,  4]
    tbus = [2, 5, 6, 7, 3, 8, 9, 10, 4, 11, 12, 13]

    # Length of each branch in km (one value per branch, same order as fbus/tbus)
    lengths = [0.6, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] .* 0.5

    # Base impedance matrix per km (Ω/km)
    Zij_per_km = ComplexF64[
        0.799+0.075im  0.353+0.004im  0.353+0.004im;
        0.353+0.004im  0.799+0.075im  0.353+0.004im;
        0.353+0.004im  0.353+0.004im  0.799+0.075im]

    # Per-branch impedance matrices: scale by length
    Zij = [lengths[k] .* Zij_per_km for k in 1:length(fbus)]

    Ysh = zeros(ComplexF64, 3, 3)

    return SystemData(nN, nL,
        P_Demand_a, P_Demand_b, P_Demand_c,
        Q_Demand_a, Q_Demand_b, Q_Demand_c,
        P_Cost_G, P_Gen_lb, P_Gen_ub, Q_Gen_lb, Q_Gen_ub, Gen_Status,
        S_Sol_Gen_ub_a, S_Sol_Gen_ub_b, S_Sol_Gen_ub_c,
        V_base, V_lb, V_ub, fbus, tbus, Zij, Ysh)
end
