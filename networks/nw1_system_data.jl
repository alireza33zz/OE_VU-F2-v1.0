# ============================================================
# networks/nw1_system_data.jl
# IEEE 13-node test feeder (modified, 3-phase unbalanced)
# Nominal voltage: 400 V line-to-line  →  400/√3 V phase
# ============================================================

"""
    NW1_V_BASE

Nominal phase voltage [V] for the 13-node network (400 V LL).
"""
const NW1_V_BASE = 400.0 / sqrt(3)

"""
    build_system_data_nw1() -> SystemData

Construct the SystemData for the 13-node, 400 V test network.
"""
function build_system_data_nw1()::SystemData
    nN = 13;  nL = 12

Demand_scaling = 1.0  # Kept at 1.0 — unbalance comes from phase asymmetry,
                           # not raw magnitude. Higher scaling risks voltage violations
                           # before VUF reaches the target range.

# Heavily loading Phase A, moderating Phase B, and stripping Phase C
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
        V_lb, V_ub, fbus, tbus, Zij, Ysh)
end

# ============================================================
# 2. ADMITTANCE MATRIX
# ============================================================

"""
    build_admittance_matrix(sys) -> Matrix{ComplexF64}

Assemble the full (3·nN × 3·nN) nodal admittance matrix Y by looping over the line list.
"""
function build_admittance_matrix(sys::SystemData)::Matrix{ComplexF64}
    n3 = 3 * sys.nN
    Y  = zeros(ComplexF64, n3, n3)

    for k in 1:sys.nL
        Yl = inv(sys.Zij[k])          # k-th branch matrix
        i = sys.fbus[k];  j = sys.tbus[k]
        ri = (3i-2):(3i);  rj = (3j-2):(3j)
        @inbounds begin
            Y[ri, ri] .+= Yl;  Y[rj, rj] .+= Yl
            Y[ri, rj] .-= Yl;  Y[rj, ri] .-= Yl
        end
    end

    return Y
end

function build_impedance_matrix(
    Y::Matrix{ComplexF64};
    slack_node::Int = 1
)::Matrix{ComplexF64}

    n3  = size(Y, 1)

    # Indices of the 3 slack bus phases to remove
    slack_cols = [3*(slack_node-1)+1, 3*(slack_node-1)+2, 3*(slack_node-1)+3]
    keep       = setdiff(1:n3, slack_cols)

    # Reduced admittance matrix — guaranteed non-singular
    Y_red = Y[keep, keep]
    Z_red = inv(Y_red)

    # Embed back into full matrix (slack rows/cols remain zero)
    Z_full = zeros(ComplexF64, n3, n3)
    Z_full[keep, keep] = Z_red

    return Z_full
end

# ============================================================
# 3. VUF CONSTRAINT CONFIGURATION
# ============================================================

"""
    VUFConfig

Specifies where and how tightly VUF is enforced.
"""
struct VUFConfig
    buses::Vector{Int}
    threshold::Float64
    pos_seq_floor::Float64

    function VUFConfig(buses::Vector{Int};
                       threshold::Float64    = 0.02,
                       pos_seq_floor::Float64 = (0.5 * 400.0/sqrt(3))^2)
        @assert 0.0 < threshold < 1.0  "VUF threshold must be in (0, 1)."
        new(buses, threshold, pos_seq_floor)
    end
end
