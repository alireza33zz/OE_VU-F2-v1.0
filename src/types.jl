# ============================================================
# types.jl — Shared struct definitions
# All power values in W / var, voltages in V, impedances in Ω.
# ============================================================

import MathOptInterface as MOI

# ── Network parameters ──────────────────────────────────────
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

    V_lb::Vector{Float64};  V_ub::Vector{Float64}

    fbus::Vector{Int};  tbus::Vector{Int}
    Zij::Vector{Matrix{ComplexF64}}   # ← changed from Matrix to Vector{Matrix}
    Ysh::Matrix{ComplexF64}
end

# ── VUF constraint configuration ────────────────────────────
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

# Convenience constructor: no VUF enforcement
VUFConfig() = VUFConfig(Int[]; threshold = 0.02)

# ── OPF solution snapshot ────────────────────────────────────
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

# ── PV unit descriptor ───────────────────────────────────────
struct PVUnit
    node::Int
    phase::Symbol
    S_max::Float64
end

# ── Operating-envelope boundary ──────────────────────────────
struct OEBoundary
    unit::PVUnit
    P_boundary::Vector{Float64}
    Q_boundary::Vector{Float64}
    P_baseline::Float64
    Q_baseline::Float64
end
