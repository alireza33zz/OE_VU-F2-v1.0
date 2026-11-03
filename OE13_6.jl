# ============================================================
# NW1 — 13-bus network
# PQ Operating Envelopes: M1 (weighted fairness) vs M2 (lexicographic max-min),
# each without and with VUF constraints.
# ============================================================

include(joinpath(@__DIR__, "src", "OEFramework.jl"))
using .OEFramework
using Plots

# ── Settings ────────────────────────────────────────────────────────────────
TEST_NUM      = 7
vuf_threshold = 0.010           # 1 % VUF limit (stricter than the 2 % IEC limit)
vuf_buses     = [11, 12, 13]    # critical buses
Angle_step    = 5               # 5° step → 72 angles
number_steps  = Int(360 / Angle_step)

# Solver settings
m2_kw = (stage_max_iter = 3000, final_max_iter = 5000)   # lex max-min Ipopt caps (default 500 / 2000)

mkpath("Graphics")
sys = build_system_data_13bus()
export_admittance_xlsx(build_admittance_matrix(sys), "system_matrix.xlsx")

# ── Warm-up calls (cheap, triggers JIT compilation) ─────────────────────────
run_case(sys; method = :weighted,  n_oe_angles = 1, show_results = false);
run_case(sys; method = :lexmaxmin, n_oe_angles = 1, show_results = false, oe_kw = m2_kw);

# ── Scenarios ───────────────────────────────────────────────────────────────
vuf_tag = "vuf_$(vuf_threshold*100)pct_TEST$(TEST_NUM)"
common  = (vuf_threshold = vuf_threshold, n_oe_angles = number_steps)

r11 = run_scenario(sys, "NW1-M1 | No VUF"; method = :weighted, common...,
                   out_dir = "opf_fairness_figures/baseline_TEST$(TEST_NUM)")
r12 = run_scenario(sys, "NW1-M1 | VUF $(vuf_threshold*100)%"; method = :weighted, common...,
                   vuf_buses = vuf_buses, out_dir = "opf_fairness_figures/$vuf_tag")
r21 = run_scenario(sys, "NW1-M2 | No VUF"; method = :lexmaxmin, oe_kw = m2_kw, common...,
                   out_dir = "opf_fairness_LEX_figures/baseline_TEST$(TEST_NUM)")
r22 = run_scenario(sys, "NW1-M2 | VUF $(vuf_threshold*100)%"; method = :lexmaxmin, oe_kw = m2_kw, common...,
                   vuf_buses = vuf_buses, out_dir = "opf_fairness_LEX_figures/$vuf_tag")

print_runtime_summary("Runtime Summary (NW1)", [
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

p_grid, p_leg = plot_operating_envelopes_overlay_legend(oe_lists)
display(p_grid)
display(p_leg)
savefig(p_grid, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_grid.svg")
savefig(p_leg,  "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_legend.svg")
