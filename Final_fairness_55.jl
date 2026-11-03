# ============================================================
# NW2 — 116-bus network
# PQ Operating Envelopes: M1 (weighted fairness) vs M2 (lexicographic max-min),
# each without and with VUF constraints.
# ============================================================

include(joinpath(@__DIR__, "src", "OEFramework.jl"))
using .OEFramework
using Plots

# ── Settings ────────────────────────────────────────────────────────────────
TEST_NUM      = 4
vuf_threshold = 0.015               # 1.5 % VUF limit (stricter than the 2 % IEC limit)
vuf_buses     = [40, 60, 80, 96]    # buses with 3-phase balanced loads
Angle_step    = 5                   # 5° step → 72 angles
number_steps  = Int(360 / Angle_step)

# Solver / plot settings specific to this (larger) network
opf_kw  = (max_iter = 10_000, tol = 1e-4)                 # baseline OPF
m1_kw   = (voltage_ref = :unit,)                          # sensitivity weights use V at the unit's bus
m2_kw   = (stage_max_iter = 3000, final_max_iter = 5000)  # lex max-min Ipopt caps
plot_kw = (nrows = 3,)

mkpath("Graphics")
sys = build_system_data_116bus()

# ── Warm-up calls (cheap, triggers JIT compilation) ─────────────────────────
run_case(sys; method = :weighted,  n_oe_angles = 1, show_results = false,
         opf_kw = opf_kw, oe_kw = m1_kw, plot_kw = plot_kw);
run_case(sys; method = :lexmaxmin, n_oe_angles = 1, show_results = false,
         opf_kw = opf_kw, oe_kw = m2_kw, plot_kw = plot_kw);

# ── Scenarios ───────────────────────────────────────────────────────────────
vuf_tag = "vuf_$(vuf_threshold*100)pct_TEST$(TEST_NUM)"
common  = (vuf_threshold = vuf_threshold, n_oe_angles = number_steps,
           opf_kw = opf_kw, plot_kw = plot_kw)

r11 = run_scenario(sys, "NW2-M1 | No VUF"; method = :weighted, oe_kw = m1_kw, common...,
                   out_dir = "opf_fairness_figures/baseline_TEST$(TEST_NUM)")
r12 = run_scenario(sys, "NW2-M1 | VUF $(vuf_threshold*100)%"; method = :weighted, oe_kw = m1_kw, common...,
                   vuf_buses = vuf_buses, out_dir = "opf_fairness_figures/$vuf_tag")
r21 = run_scenario(sys, "NW2-M2 | No VUF"; method = :lexmaxmin, oe_kw = m2_kw, common...,
                   out_dir = "opf_fairness_LEX_figures/baseline_TEST$(TEST_NUM)")
r22 = run_scenario(sys, "NW2-M2 | VUF $(vuf_threshold*100)%"; method = :lexmaxmin, oe_kw = m2_kw, common...,
                   vuf_buses = vuf_buses, out_dir = "opf_fairness_LEX_figures/$vuf_tag")

print_runtime_summary("Runtime Summary (NW2)", [
    "M1 | No VUF"                     => r11.time,
    "M1 | VUF $(vuf_threshold*100)%"  => r12.time,
    "M2 | No VUF"                     => r21.time,
    "M2 | VUF $(vuf_threshold*100)%"  => r22.time,
])

# ── Utilitarian benchmark + Price of Fairness ───────────────────────────────
# Reuses each case's baseline and M1/M2 envelopes (same F, baseline and θ_m).
u0 = run_price_of_fairness(sys, r11, r21; n_oe_angles = number_steps)
u1 = run_price_of_fairness(sys, r12, r22; vuf_buses = vuf_buses,
                           vuf_threshold = vuf_threshold, n_oe_angles = number_steps)
print_pof_table(["No VUF" => u0, "VUF $(vuf_threshold*100)%" => u1])

# ── Overlay of the four scenarios ───────────────────────────────────────────
oe_lists = [r11.oe_list, r12.oe_list, r21.oe_list, r22.oe_list]

p = plot_operating_envelopes_overlay(oe_lists)
display(p)
savefig(p, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM).pdf")

# Compact layout for the 27-unit grid (IEEE 2-column figure)
p_grid, p_leg = plot_operating_envelopes_overlay_legend(oe_lists; compact = true)
display(p_grid)
display(p_leg)
savefig(p_grid, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_grid.svg")
savefig(p_leg,  "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_legend.svg")
