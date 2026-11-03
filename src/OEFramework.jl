# ============================================================
# 3-Phase Unbalanced AC-OPF with Dynamic PQ Operating Envelopes
# and Voltage Unbalance Factor (VUF) Constraints
#
#   M1: sensitivity-weighted fairness
#   M2: lexicographic max-min fairness
#   Utilitarian benchmark → price of fairness (PoF) of M1 / M2
#
# Usage from a run script:
#   include(joinpath(@__DIR__, "src", "OEFramework.jl"))
#   using .OEFramework
# ============================================================

module OEFramework

using JuMP, Ipopt, DataFrames
import MathOptInterface as MOI
using Plots; gr()
using Statistics: mean
using LinearAlgebra, Printf, Colors
using XLSX

include("types.jl")
include("network.jl")
include("opf.jl")
include("envelopes.jl")
include("fairness.jl")
include("reporting.jl")
include("plotting.jl")
include("dashboard.jl")
include("pipeline.jl")
include("networks/nw1_13bus.jl")
include("networks/nw2_116bus.jl")

# Types
export SystemData, VUFConfig, OPFResult, PVUnit, OEBoundary
# Networks
export build_system_data_13bus, build_system_data_116bus
# Network matrices
export build_admittance_matrix, build_impedance_matrix, export_admittance_xlsx
# OPF
export build_base_model, build_and_solve_opf, extract_result, compute_vuf, compute_line_flows
# Operating envelopes
export find_pv_units, compute_sensitivity_weights_cached,
       compute_operating_envelopes, compute_operating_envelopes_lexmaxmin
# Utilitarian benchmark / price of fairness
export compute_operating_envelopes_utilitarian, aggregate_steps, polygon_area,
       price_of_fairness
# Reporting & plots
export print_results, plot_dashboard, plot_operating_envelopes,
       plot_operating_envelopes_overlay, plot_operating_envelopes_overlay_legend
# Drivers
export run_case, run_scenario, print_runtime_summary,
       run_price_of_fairness, print_pof_table

end # module
