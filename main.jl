# ============================================================
# main.jl — Run OE computation for NW1 and NW2
# ============================================================

using JuMP, Ipopt, DataFrames
import MathOptInterface as MOI
using Plots; gr()
using LinearAlgebra, Printf, Colors, XLSX
using Statistics: mean

# ── Load shared source ──────────────────────────────────────
include("src/types.jl")
include("src/admittance.jl")
include("src/opf_model.jl")
include("src/postprocess.jl")
include("src/operating_envelopes.jl")
include("src/plotting.jl")

mkpath("outputs/NW1")
mkpath("outputs/NW2")

# ── Global config ───────────────────────────────────────────
const VUF_THRESHOLD = 0.015
const N_ANGLES      = 36

# ============================================================
# Generic runner
# ============================================================
function run_case(
    name::String,
    network_file::String,
    out_dir::String;
    vuf_buses::Vector{Int}  = Int[],
    vuf_threshold::Float64  = VUF_THRESHOLD,
    vuf_V_base::Float64     = 400.0/sqrt(3),   # only used for VUFConfig floor
    n_angles::Int           = N_ANGLES
)
    println("\n", "="^60)
    println("  CASE: $name")
    println("="^60)

    include(network_file)
    sys = build_system_data()
    Y   = build_admittance_matrix(sys)

    # VUFConfig: pos_seq_floor derived from this network's V_base
    pos_seq_floor = (0.5 * vuf_V_base)^2
    vuf_cfg = VUFConfig(vuf_buses;
                        threshold     = vuf_threshold,
                        pos_seq_floor = pos_seq_floor)

    # ── Baseline OPF ──────────────────────────────────────
    println("[ 1/3 ] Baseline AC-OPF …")
    baseline_model = build_and_solve_opf(sys, Y, vuf_cfg)
    baseline_res   = extract_result(baseline_model, sys.nN)
    print_results(sys, baseline_res, vuf_buses)

    # ── Dashboard figures ──────────────────────────────────
    println("\n[ 2/3 ] Dashboard figures …")
    figs = plot_dashboard(sys, baseline_res;
                          vuf_threshold = vuf_threshold,
                          save_dir      = out_dir,
                          fmt           = "pdf")

    # ── OE — Method 1: Proportional fairness ──────────────
    println("\n[ 3/3 ] OE Method 1 (proportional fairness) …")
    t1 = @elapsed oe_m1 = compute_operating_envelopes(
    sys, Y, vuf_cfg, baseline_res; n_angles=n_angles)
    @printf "  M1 done in %.1f s\n" t1

    fig_m1 = plot_operating_envelopes(oe_m1)
    savefig(fig_m1, joinpath(out_dir, "OE_M1_subplots.pdf"))

    # ── OE — Method 2: Lex max-min ────────────────────────
    println("\n[ 3/3 ] OE Method 2 (lex max-min) …")
    t2 = @elapsed oe_m2 = compute_operating_envelopes_lexmaxmin(
    sys, Y, vuf_cfg, baseline_res; n_angles=n_angles)
    @printf "  M2 done in %.1f s\n" t2

    fig_m2 = plot_operating_envelopes(oe_m2)
    savefig(fig_m2, joinpath(out_dir, "OE_M2_subplots.pdf"))

    # ── Comparison overlay ─────────────────────────────────
    fig_cmp_grid, _ = plot_operating_envelopes_overlay_legend(
    [oe_m1, oe_m2];
    names      = ["M1 — Proportional Fairness", "M2 — Lex Max-Min"],
    colors     = [:dodgerblue, :darkorange],
    fillalphas = [0.18, 0.10],
    linestyles = [:solid, :dash]
)

    savefig(fig_cmp_grid, joinpath(out_dir, "OE_comparison.pdf"))

    println("  Figures saved → $out_dir/")
    return sys, baseline_res, oe_m1, oe_m2
end

# ============================================================
# JIT warm-up (1 angle, cheap)
# ============================================================
println("\n[ Warm-up ] JIT compilation pass …")
include("networks/nw1_system_data.jl")
_sys_w = build_system_data()
_Y_w   = build_admittance_matrix(_sys_w)
_vuf_w = VUFConfig()                         # no VUF, default floor
_m_w   = build_and_solve_opf(_sys_w, _Y_w, _vuf_w)
_r_w   = extract_result(_m_w, _sys_w.nN)
compute_operating_envelopes(_sys_w, _Y_w, _vuf_w, _r_w; n_angles=1)
compute_operating_envelopes_lexmaxmin(_sys_w, _Y_w, _vuf_w, _r_w; n_angles=1)
println("[ Warm-up done ]")

# ============================================================
# NW1 — 13-node, 400 V
# ============================================================
sys1, res1, oe1_m1, oe1_m2 = run_case(
    "NW1 — 13-node 400 V",
    "networks/nw1_system_data.jl",
    "outputs/NW1";
    vuf_buses     = NW1_VUF_BUSES,
    vuf_V_base    = NW1_V_BASE,
    n_angles      = N_ANGLES
)

include("networks/nw2_system_data.jl")

# ============================================================
# NW2 — 116-node EULV, 416 V
# ============================================================
sys2, res2, oe2_m1, oe2_m2 = run_case(
    "NW2 — 116-node EULV 416 V",
    "networks/nw2_system_data.jl",
    "outputs/NW2";
    vuf_buses     = NW2_VUF_BUSES,
    vuf_threshold = VUF_THRESHOLD,
    vuf_V_base    = NW2_V_BASE,
    n_angles      = N_ANGLES
)

println("\n✓ All cases complete. Figures in outputs/NW1/ and outputs/NW2/")