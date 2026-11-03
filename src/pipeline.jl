# ============================================================
# DRIVERS
# ============================================================

"""
    run_case(sys; method, vuf_buses, vuf_threshold, n_oe_angles,
             opf_kw, oe_kw, plot_kw, show_results) -> NamedTuple

Baseline OPF + operating envelopes for one configuration.

- `method`  : `:weighted` (M1, sensitivity-weighted fairness) or `:lexmaxmin` (M2)
- `vuf_buses`: buses where VUF ≤ `vuf_threshold` is enforced (empty = none)
- `opf_kw`  : extra keywords for `build_and_solve_opf` (e.g. `(max_iter=10_000, tol=1e-4)`)
- `oe_kw`   : extra keywords for the OE method (e.g. `(voltage_ref=:unit,)` for M1,
              `(stage_max_iter=3000, final_max_iter=5000)` for M2)
- `plot_kw` : extra keywords for `plot_operating_envelopes` (e.g. `(nrows=3,)`)

Returns `(; sys, baseline_res, oe_list, fig)`.
"""
function run_case(sys::SystemData;
    method::Symbol         = :weighted,
    vuf_buses::Vector{Int} = Int[],
    vuf_threshold::Float64 = 0.02,
    n_oe_angles::Int       = 36,
    opf_kw::NamedTuple     = (;),
    oe_kw::NamedTuple      = (;),
    plot_kw::NamedTuple    = (;),
    show_results::Bool     = true
)
    @assert method in (:weighted, :lexmaxmin)  "method must be :weighted or :lexmaxmin"

    println("[ 1/4 ]  Building admittance matrix …")
    Y = build_admittance_matrix(sys)

    vuf_cfg = VUFConfig(vuf_buses; threshold = vuf_threshold, V_nom = sys.V_nom)
    if isempty(vuf_buses)
        println("[ 2/4 ]  No VUF constraints requested.")
    else
        println("[ 2/4 ]  VUF ≤ $(vuf_threshold*100)% enforced at buses: $vuf_buses")
    end

    println("[ 3/4 ]  Solving baseline AC-OPF …")
    baseline_model = build_and_solve_opf(sys, Y, vuf_cfg; opf_kw...)
    baseline_res   = extract_result(baseline_model, sys.nN)
    show_results && print_results(sys, baseline_res, vuf_buses)

    println("\n[ 4/4 ]  Computing PQ Operating Envelopes ($method) …")
    oe_fn   = method == :weighted ? compute_operating_envelopes :
                                    compute_operating_envelopes_lexmaxmin
    oe_list = oe_fn(sys, Y, vuf_cfg, baseline_res; n_angles = n_oe_angles, oe_kw...)

    fig = plot_operating_envelopes(oe_list; plot_kw...)

    return (; sys, baseline_res, oe_list, fig)
end

"""
    run_scenario(sys, label; out_dir, vuf_threshold, fmt, kwargs...) -> NamedTuple

`run_case` (timed) + baseline dashboard + OE figure, all saved to `out_dir`.
`vuf_threshold` is also used to colour the dashboard VUF bars.
Returns the `run_case` result plus `time` [s].
"""
function run_scenario(sys::SystemData, label::String;
    out_dir::String,
    vuf_threshold::Float64 = 0.02,
    fmt::String            = "pdf",
    kwargs...
)
    println("\n", "═"^72, "\n  $label\n", "═"^72)

    t = @elapsed out = run_case(sys; vuf_threshold = vuf_threshold, kwargs...)
    t = round(t; digits = 2)
    println("► $label : $(t) s")

    println("\n[ Dashboard ]")
    figs = plot_dashboard(sys, out.baseline_res;
                          vuf_threshold = vuf_threshold,
                          save_dir      = out_dir,
                          fmt           = fmt)
    for key in (:P_power, :Q_power, :voltage, :lines)
        display(figs[key])
    end

    display(out.fig)
    savefig(out.fig, joinpath(out_dir, "pq_OE.$fmt"))

    return (; out..., time = t)
end

"""
    print_runtime_summary(title, rows)

`rows` is a vector of `label => seconds` pairs.
"""
function print_runtime_summary(title::String, rows)
    w = max(maximum(length(first(r)) for r in rows) + 2, length(title) - 14)
    line(l, r) = println("│ ", rpad(l, w), "│ ", lpad(r, 10), " s │")
    println("\n┌", "─"^(w+1), "┬", "─"^14, "┐")
    println("│ ", rpad(title, w + 15), "│")
    println("├", "─"^(w+1), "┼", "─"^14, "┤")
    for (lbl, t) in rows
        line(lbl, string(t))
    end
    println("└", "─"^(w+1), "┴", "─"^14, "┘")
end

"""
    run_price_of_fairness(sys, r_m1, r_m2; vuf_buses, vuf_threshold, n_oe_angles) -> NamedTuple

Utilitarian benchmark for one VUF case, reusing the baseline and envelopes of the
M1 / M2 results (`run_scenario` outputs) of that same case, then PoF for both.
Returns `(; oe_list, U_star, time, M1, M2)` where M1/M2 are the PoF values.
"""
function run_price_of_fairness(sys::SystemData, r_m1, r_m2;
    vuf_buses::Vector{Int} = Int[],
    vuf_threshold::Float64 = 0.02,
    n_oe_angles::Int       = 36
)
    # M1 and M2 of the same case must share the baseline (same F and starting point)
    @assert isapprox(r_m1.baseline_res.obj, r_m2.baseline_res.obj; rtol=1e-8) "M1/M2 baselines differ"

    Y   = build_admittance_matrix(sys)
    cfg = VUFConfig(vuf_buses; threshold = vuf_threshold, V_nom = sys.V_nom)

    t = @elapsed oe_util, U_star = compute_operating_envelopes_utilitarian(
            sys, Y, cfg, r_m1.baseline_res;
            n_angles = n_oe_angles, warm_lists = [r_m1.oe_list, r_m2.oe_list])

    return (; oe_list = oe_util, U_star, time = round(t; digits = 2),
              M1 = price_of_fairness(r_m1.oe_list, oe_util, U_star, n_oe_angles),
              M2 = price_of_fairness(r_m2.oe_list, oe_util, U_star, n_oe_angles))
end

"""
    print_pof_table(cases)

`cases` is a vector of `case_label => run_price_of_fairness(...)` pairs.
"""
function print_pof_table(cases)
    println("\n  Price of Fairness  (PoF = 1 − ΣA_fair / ΣA_util)")
    println("  Case       |    M1    |    M2    | util. runtime")
    for (c, u) in cases
        @printf("  %-10s | %7.2f%% | %7.2f%% | %8.2f s\n", c, 100u.M1, 100u.M2, u.time)
    end
end
