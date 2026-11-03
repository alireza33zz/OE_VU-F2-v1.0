# ============================================================
# DATA TYPES
# ============================================================

"""
    SystemData

Immutable struct holding all network parameters.
Power in W / var, voltages in V, impedances in Ω.
`V_nom` is the nominal phase-to-neutral voltage (slack reference and flat start).
"""
struct SystemData
    nN::Int
    nL::Int

    P_Demand_a::Vector{Float64};  P_Demand_b::Vector{Float64};  P_Demand_c::Vector{Float64}
    Q_Demand_a::Vector{Float64};  Q_Demand_b::Vector{Float64};  Q_Demand_c::Vector{Float64}

    P_Cost_G::Vector{Float64}
    P_Gen_lb::Vector{Float64};  P_Gen_ub::Vector{Float64}
    Q_Gen_lb::Vector{Float64};  Q_Gen_ub::Vector{Float64}
    Gen_Status::Vector{Int}

    S_Sol_Gen_ub_a::Vector{Float64}
    S_Sol_Gen_ub_b::Vector{Float64}
    S_Sol_Gen_ub_c::Vector{Float64}

    V_nom::Float64
    V_lb::Vector{Float64};  V_ub::Vector{Float64}

    fbus::Vector{Int};  tbus::Vector{Int}
    Zij::Vector{Matrix{ComplexF64}}   # one 3×3 series impedance per branch
    Ysh::Matrix{ComplexF64}
end

"""
    VUFConfig

Specifies where and how tightly VUF is enforced.
`pos_seq_floor` defaults to (0.5·V_nom)², so pass the network's `V_nom`.
"""
struct VUFConfig
    buses::Vector{Int}
    threshold::Float64
    pos_seq_floor::Float64

    function VUFConfig(buses::Vector{Int};
                       threshold::Float64     = 0.02,
                       V_nom::Float64         = 400.0/sqrt(3),
                       pos_seq_floor::Float64 = (0.5 * V_nom)^2)
        @assert 0.0 < threshold < 1.0  "VUF threshold must be in (0, 1)."
        new(buses, threshold, pos_seq_floor)
    end
end

# Convenience constructor: no VUF enforcement
VUFConfig() = VUFConfig(Int[]; threshold=0.02)

"""
    OPFResult

Extracted solution snapshot from one OPF solve.
Decouples downstream analysis from the live JuMP model.
"""
struct OPFResult
    status::MOI.TerminationStatusCode
    obj::Float64

    V_r_a::Vector{Float64};  V_i_a::Vector{Float64}
    V_r_b::Vector{Float64};  V_i_b::Vector{Float64}
    V_r_c::Vector{Float64};  V_i_c::Vector{Float64}

    P_G_a::Vector{Float64};  Q_G_a::Vector{Float64}
    P_G_b::Vector{Float64};  Q_G_b::Vector{Float64}
    P_G_c::Vector{Float64};  Q_G_c::Vector{Float64}

    P_sol_a::Vector{Float64};  Q_sol_a::Vector{Float64}
    P_sol_b::Vector{Float64};  Q_sol_b::Vector{Float64}
    P_sol_c::Vector{Float64};  Q_sol_c::Vector{Float64}
end

struct PVUnit
    node::Int
    phase::Symbol
    S_max::Float64
end

struct OEBoundary
    unit::PVUnit
    P_boundary::Vector{Float64}
    Q_boundary::Vector{Float64}
    P_baseline::Float64
    Q_baseline::Float64
end

# ── Small accessors used throughout (phase symbol → field) ──────────────────
_pick(ph::Symbol, a, b, c) = ph == :a ? a : ph == :b ? b : c

P_sol_of(res::OPFResult, u::PVUnit) = _pick(u.phase, res.P_sol_a, res.P_sol_b, res.P_sol_c)[u.node]
Q_sol_of(res::OPFResult, u::PVUnit) = _pick(u.phase, res.Q_sol_a, res.Q_sol_b, res.Q_sol_c)[u.node]
S_cap_of(sys::SystemData, u::PVUnit) =
    _pick(u.phase, sys.S_Sol_Gen_ub_a, sys.S_Sol_Gen_ub_b, sys.S_Sol_Gen_ub_c)[u.node]
