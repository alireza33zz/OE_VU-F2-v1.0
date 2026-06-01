# ============================================================
# main.jl — Top-level driver: runs NW1 and NW2, exports figures
# ============================================================
#
# Usage (REPL, from repo root):
#   include("main.jl")
#
# Terminal:
#   julia --project=. main.jl
#
# Outputs written to:
#   outputs/NW1/   <- 13-node results
#   outputs/NW2/   <- EULV 116-node results
# ============================================================

using JuMP, Ipopt, DataFrames
import MathOptInterface as MOI
using Plots; gr()
using Statistics: mean
using LinearAlgebra, Printf, Colors
using XLSX

# ── Load shared modules ──────────────────────────────────────
include("src/types.jl")
include("src/admittance.jl")
include("src/opf_model.jl")
include("src/postprocess.jl")
include("src/operating_envelopes.jl")
include("src/plotting.jl")

# ── Load network data ────────────────────────────────────────
include("networks/nw1_system_data.jl")
include("networks/nw2_system_data.jl")

# ============================================================
# Generic case runner
# ============================================================

"""
    run_case(label, sys, V_base; kwargs...) -> NamedTuple

Run a full OE study for one test case:
  1. Build admittance matrix
  2. Solve baseline AC-OPF
  3. Compute OEs — weighted fairness and lexicographic max-min
  4. Export all figures to `out_dir`

Returns a NamedTuple: sys, Y, baseline_res, oe_fair, oe_lex,
                      figs_dash, fig_oe_fair, fig_oe_lex
"""
function run_case(
    label::String,
    sys::SystemData,
    V_base::Float64;
    vuf_buses::Vector{Int} = Int[],
    vuf_threshold::Float64 = 0.02,
    n_oe_angles::Int       = 36,
    out_dir::String        = joinpath("outputs", label)
)
    mkpath(out_dir)
    println("\n", "="^60)
    println("  Case: $label")
    println("="^60)

    # Step 1 — admittance matrix
    println("[ 1/5 ] Building admittance matrix ...")
    Y = build_admittance_matrix(sys)

    # Step 2 — VUF configuration
    vuf_cfg = VUFConfig(vuf_buses; threshold = vuf_threshold)
    if isempty(vuf_buses)
        println("[ 2/5 ] No VUF constraints.")
    else
        println("[ 2/5 ] VUF <= $(vuf_threshold*100)% at buses: $vuf_buses")
    end

    # Step 3 — baseline AC-OPF
    println("[ 3/5 ] Solving baseline AC-OPF ...")
    bm           = build_and_solve_opf(sys, Y, vuf_cfg; V_base = V_base)
    baseline_res = extract_result(bm, sys.nN)
    print_results(sys, baseline_res, vuf_buses)

    # Step 4 — operating envelopes
    println("[ 4/5 ] Computing OEs (weighted fairness) ...")
    oe_fair = compute_operating_envelopes(
        sys, Y, vuf_cfg, baseline_res;
        n_angles = n_oe_angles,
        V_base   = V_base)

    println("[ 4/5 ] Computing OEs (lexicographic max-min) ...")
    oe_lex = compute_operating_envelopes_lexmaxmin(
        sys, Y, vuf_cfg, baseline_res;
        n_angles = n_oe_angles,
        V_base   = V_base)

    # Step 5 — figures
    println("[ 5/5 ] Generating figures -> $out_dir")

    figs_dash = plot_dashboard(sys, baseline_res;
        vuf_threshold = vuf_threshold,
        save_dir      = out_dir,
        fmt           = "pdf")

    fig_oe_fair = plot_operating_envelopes_overlay_legend(oe_fair)
    savefig(fig_oe_fair, joinpath(out_dir, "oe_fairness.pdf"))
    println("  Saved -> $(joinpath(out_dir, "oe_fairness.pdf"))")

    fig_oe_lex = plot_operating_envelopes_overlay_legend(oe_lex)
    savefig(fig_oe_lex, joinpath(out_dir, "oe_lexmaxmin.pdf"))
    println("  Saved -> $(joinpath(out_dir, "oe_lexmaxmin.pdf"))")

    if length(oe_fair) == length(oe_lex)
        fig_cmp = plot_operating_envelopes_overlay(
            [oe_fair, oe_lex], ["Fairness", "LexMaxMin"])
        savefig(fig_cmp, joinpath(out_dir, "oe_comparison.pdf"))
        println("  Saved -> $(joinpath(out_dir, "oe_comparison.pdf"))")
    end

    println("  Case $label complete.")
    return (sys          = sys,
            Y            = Y,
            baseline_res = baseline_res,
            oe_fair      = oe_fair,
            oe_lex       = oe_lex,
            figs_dash    = figs_dash,
            fig_oe_fair  = fig_oe_fair,
            fig_oe_lex   = fig_oe_lex)
end

# ============================================================
# JIT warm-up  (cheap single-angle run to compile all methods)
# ============================================================
println("\n[ JIT warm-up ] ...")
let
    _sys = build_system_data_nw1()
    _Y   = build_admittance_matrix(_sys)
    _cfg = VUFConfig()
    _bm  = build_and_solve_opf(_sys, _Y, _cfg;
               V_base   = NW1_V_BASE,
               max_iter = 500,
               tol      = 1e-3)
    _res = extract_result(_bm, _sys.nN)
    compute_operating_envelopes(_sys, _Y, _cfg, _res;
        n_angles = 1, V_base = NW1_V_BASE)
    compute_operating_envelopes_lexmaxmin(_sys, _Y, _cfg, _res;
        n_angles = 1, V_base = NW1_V_BASE)
end
println("[ JIT warm-up ] done.\n")

# ============================================================
# Study parameters  — edit these before running
# ============================================================
VUF_THRESHOLD = 0.010   # 1.0 % VUF limit (stricter than 2 % IEC default)
ANGLE_STEP    = 5       # degrees per tracing step  =>  72 directions
N_ANGLES      = 360 ÷ ANGLE_STEP

# ============================================================
# NW1 — IEEE 13-node, 400 V
# ============================================================
t_nw1 = @elapsed begin
    nw1_result = run_case(
        "NW1",
        build_system_data_nw1(),
        NW1_V_BASE;
        vuf_buses     = [11, 12, 13],
        vuf_threshold = VUF_THRESHOLD,
        n_oe_angles   = N_ANGLES
    )
end
@printf("\n NW1 total time : %.2f s\n", t_nw1)

# ============================================================
# NW2 — EULV IEEE 116-node, 416 V
# ============================================================
t_nw2 = @elapsed begin
    nw2_result = run_case(
        "NW2",
        build_system_data_nw2(),
        NW2_V_BASE;
        vuf_buses     = collect(80:116),
        vuf_threshold = VUF_THRESHOLD,
        n_oe_angles   = N_ANGLES
    )
end
@printf("\n NW2 total time : %.2f s\n", t_nw2)

println("\nAll cases complete.  Results in outputs/")
