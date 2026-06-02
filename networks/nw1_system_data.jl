"""
build_system_data() -> SystemData
13-node, 400 V radial test network (NW1).
"""
function build_system_data()::SystemData
    nN = 13; nL = 12

    P_Demand_a = zeros(nN); P_Demand_b = zeros(nN); P_Demand_c = zeros(nN)
    Q_Demand_a = zeros(nN); Q_Demand_b = zeros(nN); Q_Demand_c = zeros(nN)

    for b in 3:13
        P_Demand_a[b] = 3000.0; P_Demand_b[b] = 2000.0; P_Demand_c[b] = 1000.0
    end
    Q_Demand_a = P_Demand_a .* 0.25
    Q_Demand_b = P_Demand_b .* 0.25
    Q_Demand_c = P_Demand_c .* 0.25

    S_Sol_Gen_ub_a = zeros(nN); S_Sol_Gen_ub_b = zeros(nN); S_Sol_Gen_ub_c = zeros(nN)
    for b in [3, 5, 7, 9, 11, 13]
        S_Sol_Gen_ub_a[b] = 7500.0; S_Sol_Gen_ub_b[b] = 7500.0; S_Sol_Gen_ub_c[b] = 7500.0
    end

    P_Cost_G   = [1.0; zeros(nN-1)]
    P_Gen_lb   = [-5_000_000.0; zeros(nN-1)]; P_Gen_ub = [5_000_000.0; zeros(nN-1)]
    Q_Gen_lb   = [-1_000_000.0; zeros(nN-1)]; Q_Gen_ub = [1_000_000.0; zeros(nN-1)]
    Gen_Status = [1; zeros(Int, nN-1)]

    V_base = 400.0 / sqrt(3)
    V_lb   = fill(0.94 * V_base, nN); V_ub = fill(1.10 * V_base, nN)

    fbus = [1, 2, 2, 3,  4,  5,  6,  7,  8,  9, 10, 11]
    tbus = [2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13]

    r1 = 0.411; x1 = 0.272; r0 = 1.107; x0 = 0.272
    z1 = r1+im*x1; z0 = r0+im*x0
    zs = (z0+2*z1)/3; zm = (z0-z1)/3
    Z_per_km = [zs zm zm; zm zs zm; zm zm zs]
    Zij = [Z_per_km .* 0.035 for _ in 1:nL]   # 35 m branches
    Ysh = zeros(ComplexF64, 3, 3)

    return SystemData(nN, nL,
        P_Demand_a, P_Demand_b, P_Demand_c,
        Q_Demand_a, Q_Demand_b, Q_Demand_c,
        P_Cost_G, P_Gen_lb, P_Gen_ub, Q_Gen_lb, Q_Gen_ub, Gen_Status,
        S_Sol_Gen_ub_a, S_Sol_Gen_ub_b, S_Sol_Gen_ub_c,
        V_lb, V_ub, fbus, tbus, Zij, Ysh)
end

const NW1_V_BASE    = 400.0 / sqrt(3)
const NW1_VUF_BUSES = [11, 12, 13]
const NW1_NAME      = "NW1_13nodes"