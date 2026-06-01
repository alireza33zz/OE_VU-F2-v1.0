# ============================================================
# SHARED DATA STRUCTURES
# ============================================================
# This module defines core data structures used across all network analysis scenarios.
# All structures are immutable for thread-safety and clarity of intent.

using JuMP
import MathOptInterface as MOI

"""
    SystemData

Immutable struct holding all network parameters.
Power in W / var, voltages in V, impedances in Ω.

Fields:
- nN::Int              : Number of nodes
- nL::Int              : Number of lines
- P_Demand_*::Vector   : Real power demand per phase [W]
- Q_Demand_*::Vector   : Reactive power demand per phase [var]
- P_Cost_G::Vector     : Generation cost coefficients [£/W] or similar
- P_Gen_*::Vector      : Real power generation bounds [W]
- Q_Gen_*::Vector      : Reactive power generation bounds [var]
- Gen_Status::Vector   : Commitment status (0=offline, 1=online)
- S_Sol_Gen_ub_*::Vector : Maximum PV capacity per phase [VA]
- V_*::Vector          : Voltage bounds [V] (line-to-neutral)
- fbus, tbus::Vector   : From/to bus indices for each line
- Zij::Vector{Matrix}  : Per-branch impedance matrix (3×3 complex) [Ω]
- Ysh::Matrix          : Shunt admittance matrix (typically sparse/zero)
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

    V_lb::Vector{Float64};  V_ub::Vector{Float64}

    fbus::Vector{Int};  tbus::Vector{Int}
    Zij::Vector{Matrix{ComplexF64}}
    Ysh::Matrix{ComplexF64}
end

"""
    VUFConfig

Configuration for Voltage Unbalance Factor (VUF) constraints.
Specifies where and how tightly VUF is enforced during optimization.

Fields:
- buses::Vector{Int}              : Bus indices where VUF is constrained
- threshold::Float64              : Maximum acceptable VUF (0 < ε < 1, typically 0.02 = 2%)
- pos_seq_floor::Float64          : Minimum positive-sequence voltage squared [V²]
                                   Prevents trivial solutions with |V+| → 0
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

"""
    VUFConfig()

Convenience constructor: no VUF enforcement (empty bus list).
"""
VUFConfig() = VUFConfig(Int[]; threshold=0.02)

"""
    OPFResult

Extracted solution snapshot from one OPF solve.
Decouples downstream analysis from the live JuMP model to enable
post-processing without keeping the solver object in memory.

Fields:
- status::MOI.TerminationStatusCode : Solver termination status (OPTIMAL, INFEASIBLE, etc.)
- obj::Float64                      : Objective value (cost)
- V_r_*::Vector, V_i_*::Vector     : Voltage rectangular coordinates per phase [V]
- P_G_*, Q_G_*::Vector             : Generation per phase [W, var]
- P_sol_*, Q_sol_*::Vector         : Solar injection per phase [W, var]
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

"""
    PVUnit

Identifies a solar PV generation unit at a specific bus and phase.

Fields:
- bus::Int    : Bus index where PV is located
- phase::Int  : Phase (1=A, 2=B, 3=C)
"""
struct PVUnit
    bus::Int
    phase::Int
end

"""
    OEBoundary

One boundary point of the operating envelope for a PV unit.

Fields:
- P_PV::Float64  : Real power injection [W]
- Q_PV::Float64  : Reactive power injection [var]
- θ::Float64     : Power angle (direction) used to find this point [rad]
"""
struct OEBoundary
    P_PV::Float64
    Q_PV::Float64
    θ::Float64
end
