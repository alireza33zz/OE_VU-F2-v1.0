# ============================================================
# plotting.jl — All visualisation helpers
# ============================================================

"""
    print_results(sys, res, vuf_buses)

Print all result tables to stdout.
"""
function print_results(sys::SystemData, res::OPFResult, vuf_buses::Vector{Int}=Int[])
    println("\n", "="^72)
    println("  3-PHASE AC-OPF RESULTS")
    println("="^72)
    println("  Status    : ", res.status)
    println("  Objective : ", round(res.obj; digits=2), "  (total generation cost)")
    println("="^72, "\n")

    nN  = sys.nN

    # ---- Voltage magnitudes & angles ----
    mag(r,i) = sqrt.(r.^2 .+ i.^2)
    arg(r,i) = rad2deg.(atan.(i, r))

    Va_mag = mag(res.V_r_a, res.V_i_a);  Va_arg = arg(res.V_r_a, res.V_i_a)
    Vb_mag = mag(res.V_r_b, res.V_i_b);  Vb_arg = arg(res.V_r_b, res.V_i_b)
    Vc_mag = mag(res.V_r_c, res.V_i_c);  Vc_arg = arg(res.V_r_c, res.V_i_c)

    # ---- VUF (post-processing) ----
    vuf = compute_vuf(res, nN)

    # Flag buses where VUF is constrained
    flag = [i in vuf_buses ? "*" : " " for i in 1:nN]

    println("── Nodal voltages — magnitude (V), angle (°), VUF (%)  [* = VUF constrained] ──")
    println(DataFrame(
        Node    = 1:nN,
        Flag    = flag,
        Va_mag  = round.(Va_mag; digits=2),
        Va_arg  = round.(Va_arg; digits=1),
        Vb_mag  = round.(Vb_mag; digits=2),
        Vb_arg  = round.(Vb_arg; digits=1),
        Vc_mag  = round.(Vc_mag; digits=2),
        Vc_arg  = round.(Vc_arg; digits=1),
        VUF_pct = round.(vuf;    digits=2)
    ))

    println("\n── Conventional generator dispatch (W / var) ──")
    println(DataFrame(
        Node = 1:nN,
        PG_a = round.(res.P_G_a; digits=2),  QG_a = round.(res.Q_G_a; digits=2),
        PG_b = round.(res.P_G_b; digits=2),  QG_b = round.(res.Q_G_b; digits=2),
        PG_c = round.(res.P_G_c; digits=2),  QG_c = round.(res.Q_G_c; digits=2)))

    println("\n── PV generator dispatch (W / var) ──")
    println(DataFrame(
        Node   = 1:nN,
        Psol_a = round.(res.P_sol_a; digits=2),  Qsol_a = round.(res.Q_sol_a; digits=2),
        Psol_b = round.(res.P_sol_b; digits=2),  Qsol_b = round.(res.Q_sol_b; digits=2),
        Psol_c = round.(res.P_sol_c; digits=2),  Qsol_c = round.(res.Q_sol_c; digits=2)))

    Pij, Pji, Qij, Qji = compute_line_flows(sys, res)

    println("\n── Active power flows (W)  [ij = sending end, ji = receiving end] ──")

"""
    plot_operating_envelopes(oe_list) -> Plots.Plot
"""
function plot_operating_envelopes(oe_list::Vector{OEBoundary})
    n     = length(oe_list)
    nrows = floor(Int, sqrt(n * 1.0))
    ncols = ceil(Int, n / nrows)
    scale = 1e-3   # W → kW

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=900)

    θ_circ = range(0, 2π, length=200)

    for (k, oe) in enumerate(oe_list)
        ax    = p[k]
        S_max = oe.unit.S_max * scale

        # Hardware S_max circle
        plot!(ax, S_max.*cos.(θ_circ), S_max.*sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label="S_max")

        # OE shaded polygon
        if length(oe.P_boundary) >= 3
            cx = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
            idx = sortperm(atan.(oe.Q_boundary.-cy, oe.P_boundary.-cx))
            Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
            Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale
            plot!(ax, Ps, Qs; seriestype=:shape,
                  fillcolor=:dodgerblue, fillalpha=0.18,
                  lc=:dodgerblue, lw=1.8, label="OE")
        end

        # Baseline operating point
        scatter!(ax, [oe.P_baseline*scale], [oe.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label="Baseline")

        # Zero-crossing lines
        hline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)
        vline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)

        title!(ax, "Node $(oe.unit.node) · Ph $(uppercase(string(oe.unit.phase)))")
        xlabel!(ax, "P (kW)");  ylabel!(ax, "Q (kvar)")
        plot!(ax, legend=(k==1), legendfont=font(7), titlefont=font(9,:bold),

# ── Dashboard figures ───────────────────────────────────────
"""
    _fig_active_power(sys, res) -> Plots.Plot

Three-column cluster per node (Ph-A, Ph-B, Ph-C). Within each phase, two bars:
  ▼  Active demand  (negative → downward)
  ▲  PV injection   (positive → upward)

x-axis: node index.  y-axis: kW.
Legend is arranged as 3 columns × 2 rows and placed outside (top center).
Solar injection uses a distinct hatch fill.
"""
function _fig_active_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3   # W → kW

    # --- build bar matrices (nN × 3) -----------------------------------------
    Pd = hcat(sys.P_Demand_a, sys.P_Demand_b, sys.P_Demand_c) .* scale
    Pp = hcat(res.P_sol_a,    res.P_sol_b,    res.P_sol_c)    .* scale

    # Clamp near-zero to exactly 0 for clean rendering
    Pd[abs.(Pd) .< 1e-2] .= 0.0
    Pp[abs.(Pp) .< 1e-2] .= 0.0

    cols_solid = [_PA, _PB, _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Demand Ph-A", "Demand Ph-B", "Demand Ph-C"]
    labels_pv  = ["PV Ph-A",     "PV Ph-B",     "PV Ph-C"]

    # Phase cluster offsets (A,B,C) → (−δ, 0, +δ)
    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)

    # Inner separation within a phase (Demand vs PV)
    δ_inner   = 0.09
    w_bar     = 0.16

    p = plot(; _FIG_KW(1320, 520)...,
        title        = "Active Power — Demand (↓) and PV Injection (↑)",
        titlefont    = (11, _FG, :bold),
        xlabel       = "Node",
        ylabel       = "P  [kW]",
        xticks       = (nodes, string.(nodes)),
        xrotation    = 0,
        ylims        = (-maximum(abs.(Pd))*1.1, maximum(Pp)*1.1),
        # Legend requested outside, top center

"""
    _fig_reactive_power(sys, res) -> Plots.Plot

Same layout as Fig 1, but for reactive power (kvar).
Three-column cluster per node; legend outside (top center), 3×2 layout.
Solar injection uses a distinct hatch fill.
"""
function _fig_reactive_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3

    Qd = hcat(sys.Q_Demand_a, sys.Q_Demand_b, sys.Q_Demand_c) .* scale
    Qp = hcat(res.Q_sol_a,    res.Q_sol_b,    res.Q_sol_c)    .* scale

    Qd[abs.(Qd) .< 1e-4] .= 0.0
    Qp[abs.(Qp) .< 1e-4] .= 0.0

    cols_solid = [_PA, _PB, _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Q-demand Ph-A", "Q-demand Ph-B", "Q-demand Ph-C"]
    labels_pv  = ["Q-PV Ph-A",     "Q-PV Ph-B",     "Q-PV Ph-C"]

    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    δ_inner   = 0.09
    w_bar     = 0.16

    p = plot(; _FIG_KW(1320, 520)...,
        title        = "Reactive Power — Demand (↓) and PV Injection (↑)",
        titlefont    = (11, _FG, :bold),
        xlabel       = "Node",
        ylabel       = "Q  [kvar]",
        xticks       = (nodes, string.(nodes)),
        ylims        = (-maximum(abs.(Qd))*1.1, maximum(Qp)*1.1),
        # Legend requested outside, top center

"""
    _fig_voltage(sys, res, vuf_threshold) -> Plots.Plot

Two vertically stacked subplots sharing the node x-axis:

  Top    : Phase voltage magnitudes |Va|, |Vb|, |Vc| [V] with ±3% limit band.
  Bottom : VUF [%] bar per node, coloured by severity.
           Nodes breaching `vuf_threshold` are drawn in _VUF_HI orange.
"""
function _fig_voltage(sys::SystemData, res::OPFResult,
                      vuf_threshold::Float64 = 0.02)

    nN    = sys.nN
    nodes = collect(_NODES)
    vuf   = compute_vuf(res, nN)

    _mag(r_, i_) = sqrt.(r_.^2 .+ i_.^2)
    Va = _mag(res.V_r_a, res.V_i_a)
    Vb = _mag(res.V_r_b, res.V_i_b)
    Vc = _mag(res.V_r_c, res.V_i_c)

    Vmin = sys.V_lb[1];  Vmax = sys.V_ub[1]
    Vnom = 400.0 / sqrt(3)

    # ── Sub-layout: 2 rows, height ratio 3:1.5 ───────────────────────────────
    p = plot(layout  = grid(2, 1; heights=[0.62, 0.38]),
             ; _FIG_KW(1320, 640)...,
             link    = :x,
             bottom_margin = 4Plots.mm,
             left_margin   = 6Plots.mm)

    # ── TOP: voltage magnitudes ───────────────────────────────────────────────
    # Limit band (filled region between Vmin and Vmax)
    plot!(p[1],
        [nodes[1]-0.5; nodes; nodes[end]+0.5;

"""
    _fig_line_flows(sys, res) -> Plots.Plot

Two vertically stacked subplots:

  Top    : Apparent power |S_ij| = √(P²+Q²) per line per phase [kVA].
           Rendered as a 3-column cluster per line (Ph-A, Ph-B, Ph-C).
           Bars are hatched when active power direction is negative (P_ij < 0).

  Bottom : Active power losses per line [W], non-negative for display:
           Sum over phases of | |P_ij| − |P_ji| |.

x-axis labels show the line as "i→j".
"""
function _fig_line_flows(sys::SystemData, res::OPFResult)
    nL = sys.nL
    Pij, Pji, Qij, Qji = compute_line_flows(sys, res)

    # Apparent power at sending end  [kVA]  (nL × 3)

"""
    plot_dashboard(sys, res;
                  vuf_threshold, save_dir, fmt) -> OrderedDict

Build and return all four baseline-OPF figures.

Figures
───────
  :P_power   — Active power demand + PV injection per node/phase
  :Q_power   — Reactive power demand + PV injection per node/phase
  :voltage   — Voltage magnitudes + VUF per node
  :lines     — Apparent power flows + losses per line

Parameters
──────────
- vuf_threshold : fraction (e.g. 0.02 = 2 %) used to colour VUF bars
- save_dir      : directory to write PNG files; "" = do not save
- fmt           : file extension (default "png")

Returns a Dict so callers can grab individual figures by key.
"""
function plot_dashboard(
        sys::SystemData,
        res::OPFResult;
        vuf_threshold::Float64 = 0.02,
        save_dir::String       = "",
        fmt::String            = "png"
    )::Dict{Symbol, Plots.Plot}

    println("  Building Fig 1/4 — Active power …")
    f1 = _fig_active_power(sys, res)

    println("  Building Fig 2/4 — Reactive power …")
    f2 = _fig_reactive_power(sys, res)

    println("  Building Fig 3/4 — Voltages + VUF …")
    f3 = _fig_voltage(sys, res, vuf_threshold)

    println("  Building Fig 4/4 — Line flows + losses …")
    f4 = _fig_line_flows(sys, res)

    figs = Dict(:P_power => f1, :Q_power => f2, :voltage => f3, :lines => f4)

    if !isempty(save_dir)
        isdir(save_dir) || mkpath(save_dir)
        for (key, fig) in figs
            path = joinpath(save_dir, "opf_$(key).$fmt")
            savefig(fig, path)
            println("  Saved → $path")
        end
    end

    return figs
end

# ── Overlay / legend figures (IEEE-formatted) ───────────────
"""
    plot_operating_envelopes_overlay(oe_lists; names, colors, fillalphas, linestyles) -> Plot

Overlays multiple OE scenarios (each is Vector{OEBoundary}) on the same grid layout.
Alignment is done by (node, phase) using oe.unit.node and oe.unit.phase.
"""
function plot_operating_envelopes_overlay(
    oe_lists::Vector{Vector{OEBoundary}};
    names      = ["Scen1","Scen2","Scen3","Scen4"],
    colors     = [:dodgerblue, :dodgerblue, :darkorange, :darkorange],
    fillalphas = [0.10, 0.18, 0.10, 0.18],
    linestyles = [:solid, :dash, :solid, :dash],
    scale      = 1e-3,
    dpi        = 900
)

    @assert length(oe_lists) == length(names) == length(colors) ==
            length(fillalphas) == length(linestyles)

    key_of(oe) = (oe.unit.node, oe.unit.phase)

    #  rename from "keys" to "ref_keys"
    ref_keys = [key_of(oe) for oe in oe_lists[1]]
    n = length(ref_keys)

    # Lookup map for each scenario
    maps = [Dict(key_of(oe) => oe for oe in oelist) for oelist in oe_lists]

    # sanity check (now uses Base.keys)
    for (i, mp) in enumerate(maps)
        missing = setdiff(ref_keys, collect(Base.keys(mp)))
        if !isempty(missing)
            error("Scenario $(i) is missing node/phase keys: $(missing)")
        end
    end

    # Layout
    nrows = 3
    ncols = ceil(Int, n / nrows)

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=dpi)

    θ_circ = range(0, 2π, length=200)

    for (k, key) in enumerate(ref_keys)
        ax = p[k]

        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        # Smax circle
        plot!(ax, S_max .* cos.(θ_circ), S_max .* sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label=false)

        # Overlay scenarios
        for s in eachindex(oe_lists)
            oe = maps[s][key]

            if length(oe.P_boundary) >= 3
                cx = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
                idx = sortperm(atan.(oe.Q_boundary .- cy, oe.P_boundary .- cx))
                Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
                Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale

                plot!(ax, Ps, Qs;
                      seriestype=:shape,
                      fillcolor=colors[s], fillalpha=fillalphas[s],
                      lc=colors[s], lw=1.8, ls=linestyles[s],
                      label = (k==1 ? names[s] : false))
            end
        end

        # Baseline star (reference scenario)
        scatter!(ax, [oe_ref.P_baseline*scale], [oe_ref.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label=false)

        hline!(ax, [0.0]; lc=:black, lw=0.6, ls=:dot, label=false)
        vline!(ax, [0.0]; lc=:black, lw=0.6, ls=:dot, label=false)

        title!(ax, "Node $(oe_ref.unit.node) · Ph $(uppercase(string(oe_ref.unit.phase)))")
        xlabel!(ax, "P (kW)"); ylabel!(ax, "Q (kvar)")

        plot!(ax, legend=(k==1), legendfont=font(7), titlefont=font(9,:bold),
              tickfont=font(7), guidefont=font(8),
              aspect_ratio=:equal, framestyle=:box, grid=true, gridalpha=0.25)
    end

    for k in n+1:nrows*ncols
        plot!(p[k], axis=false, ticks=false, grid=false, label=false)
    end

    return p
end

p = plot_operating_envelopes_overlay(
    [oe_list11, oe_list12, oe_list21, oe_list22];
    names      = ["M1 baseline", "M1 +VUF", "M2 baseline", "M2 +VUF"],
    colors     = [colorant"#56B4E9", colorant"#0072B2",
                  colorant"#E69F00", colorant"#D55E00"],
    linestyles = [:solid, :dash, :solid, :dash],
    fillalphas = [0.2, 0.05, 0.2, 0.05]
)

display(p)

savefig(p, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM).pdf")

using FileIO, ImageIO, ImageCore

font_name = "Courier"

function plot_operating_envelopes_overlay_legend(
    oe_lists::Vector{Vector{OEBoundary}};
    names      = ["Scen1","Scen2","Scen3","Scen4"],
    colors     = [:dodgerblue, :dodgerblue, :darkorange, :darkorange],
    fillalphas = [0.10, 0.18, 0.10, 0.18],
    linestyles = [:solid, :dash, :solid, :dash],
    scale      = 1e-3,
    dpi        = 900
)
    @assert length(oe_lists) == length(names) == length(colors) ==
            length(fillalphas) == length(linestyles)

    key_of(oe) = (oe.unit.node, oe.unit.phase)
    ref_keys   = [key_of(oe) for oe in oe_lists[1]]
    n          = length(ref_keys)
    maps       = [Dict(key_of(oe) => oe for oe in oelist) for oelist in oe_lists]

    for (i, mp) in enumerate(maps)
        missing_keys = setdiff(ref_keys, collect(Base.keys(mp)))
        isempty(missing_keys) || error("Scenario $i missing keys: $missing_keys")
    end

    nrows = 3
    ncols = ceil(Int, n / nrows)   # = 9 for 27 units

    # ── Wider per-cell size + larger margins for 2-col IEEE figure ───────────
    cell_w = 260
    cell_h = 260
    W = ncols * cell_w
    H = nrows * cell_h

    p_grid = plot(layout = (nrows, ncols),
                  size   = (W, H),
                  margin = 0Plots.mm,          # ← zero global margin; we control per-cell
                  background_color = :white,
                  dpi = dpi,
                  fontfamily = font_name)

    θ_circ = range(0, 2π, length=200)

    for (k, key) in enumerate(ref_keys)
        ax     = p_grid[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(θ_circ), S_max .* sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=0.8, label=false)   # thinner circle

        for s in eachindex(oe_lists)
            oe = maps[s][key]
            if length(oe.P_boundary) >= 3
                cx  = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
                idx = sortperm(atan.(oe.Q_boundary .- cy, oe.P_boundary .- cx))
                Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
                Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale

                plot!(ax, Ps, Qs;
                      seriestype = :shape,
                      fillcolor  = colors[s], fillalpha = fillalphas[s],
                      lc = colors[s], lw = 1.2, ls = linestyles[s],  # thinner lw
                      label = false)
            end
        end

        # scatter!(ax, [oe_ref.P_baseline*scale], [oe_ref.Q_baseline*scale];
        #         marker=:star5, ms=5, mc=:crimson, msc=:crimson, label=false)  # smaller star

        hline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)
        vline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)

        # ── Compact title: drop "Node"/"Ph" text to save horizontal space ───
        title!(ax, "N$(oe_ref.unit.node)·$(uppercase(string(oe_ref.unit.phase)))")

        # ── Axis labels: edge cells only, but legible fonts ─────────────────
        row = div(k - 1, ncols) + 1
        col = mod(k - 1, ncols) + 1

        xlabel!(ax, row == nrows ? "P (kW)"   : "")
        ylabel!(ax, col == 1     ? "Q (kvar)" : "")

        if col != 1
            plot!(ax, yformatter = _ -> "", ytickfontsize = 0)
        end
        if row != nrows
            plot!(ax, xformatter = _ -> "", xtickfontsize = 0)
        end

        # Asymmetric margins per cell position
        _lm = col == 1     ? 11Plots.mm  : 1Plots.mm
        _rm = 4Plots.mm    # ← increase for ALL cells, not just last col
        _bm = row == nrows ? 9Plots.mm  : 2.5Plots.mm   # ← 5mm even for interior rows
        _tm = row == 1     ? 3Plots.mm  : 2.5Plots.mm   # ← 6mm for rows 2,3 to absorb xtick bleed

        plot!(ax,
              legend         = false,
              titlefont      = font(12, font_name, :bold),
              tickfont       = font(10, font_name),
              guidefont      = font(11, font_name),
              aspect_ratio   = :equal,
              framestyle     = :box,
              grid           = true,
              gridalpha      = 0.25,
              tick_direction = :in,             # ← :in keeps ticks inside, no external bleed
              ytick_direction = :in,
              yformatter     = col == 1 ? (v -> string(Int(v))) : (_ -> ""),
              xticks         = [-6, -3, 0, 3, 6],
              yticks         = [-6, -3, 0, 3, 6],
              yflip          = false,
              ymirror        = false,           # ← stops right-side y-tick duplication
              xmirror        = false,
              left_margin    = _lm,
              right_margin   = _rm,
              top_margin     = _tm,
              bottom_margin  = _bm)
    end

    for k in n+1 : nrows*ncols
        plot!(p_grid[k], axis=false, ticks=false, grid=false, label=false)
    end

    # ── legend strip (unchanged logic, just full width) ──────────────────────
    p_leg = plot(size             = (W, 55),
                 margin           = 0Plots.mm,
                 background_color = :white,
                 dpi              = dpi,
                 fontfamily       = font_name)

    for s in eachindex(oe_lists)
        plot!(p_leg, Float64[], Float64[];
              seriestype = :shape,
              fillcolor  = colors[s], fillalpha = fillalphas[s],
              lc = colors[s], lw = 1.8, ls = linestyles[s],
              label = names[s])
    end

    plot!(p_leg,
          axis                    = false,
          ticks                   = false,
          grid                    = false,
          legend                  = :inside,
          legendfont              = font(10, font_name),
          legend_columns          = length(oe_lists),
          background_color_inside = :white,
          foreground_color_legend = :black)

    return p_grid, p_leg
end

# helper
fill_img(c, h, w) = fill(c, h, w)


p_grid, p_leg = plot_operating_envelopes_overlay_legend(
                [oe_list11, oe_list12, oe_list21, oe_list22];
                names      = ["M1 baseline", "M1 +VUF", "M2 baseline", "M2 +VUF"],
                colors     = [colorant"#56B4E9", colorant"#0072B2",
                colorant"#E69F00", colorant"#D55E00"],
                linestyles = [:solid, :dash, :solid, :dash],
                fillalphas = [0.2, 0.05, 0.2, 0.05]
)

# display in REPL / VSCode
display(p_grid)
display(p_leg)

# save — two options:

# Option A: save as PNG (raster, fine for submission if dpi=900)
savefig(p_grid, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_grid.svg")
savefig(p_leg,  "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_legend.svg")

# Option B: true single stacked PDF via pdfunite (if you have poppler installed)
# run(`pdfunite Graphics/overlay_M1_M2_states_grid.pdf Graphics/overlay_M1_M2_states_legend.pdf Graphics/overlay_M1_M2_states.pdf`)

function plot_operating_envelopes_overlay_legend(
    oe_lists::Vector{Vector{OEBoundary}};
    names      = ["Scen1","Scen2","Scen3","Scen4"],
    colors     = [:dodgerblue, :dodgerblue, :darkorange, :darkorange],
    fillalphas = [0.10, 0.18, 0.10, 0.18],
    linestyles = [:solid, :dash, :solid, :dash],
    scale      = 1e-3,
    dpi        = 900
)
    @assert length(oe_lists) == length(names) == length(colors) ==
            length(fillalphas) == length(linestyles)

    key_of(oe) = (oe.unit.node, oe.unit.phase)
    ref_keys   = [key_of(oe) for oe in oe_lists[1]]
    n          = length(ref_keys)
    maps       = [Dict(key_of(oe) => oe for oe in oelist) for oelist in oe_lists]

    for (i, mp) in enumerate(maps)
        missing_keys = setdiff(ref_keys, collect(Base.keys(mp)))
        isempty(missing_keys) || error("Scenario $i missing keys: $missing_keys")
    end

    nrows = 3
    ncols = ceil(Int, n / nrows)   # = 9 for 27 units

    # ── Wider per-cell size + larger margins for 2-col IEEE figure ───────────
    cell_w = 260
    cell_h = 260
    W = ncols * cell_w
    H = nrows * cell_h

    p_grid = plot(layout = (nrows, ncols),
                  size   = (W, H),
                  margin = 0Plots.mm,          # ← zero global margin; we control per-cell
                  background_color = :white,
                  dpi = dpi,
                  fontfamily = font_name)

    θ_circ = range(0, 2π, length=200)

    for (k, key) in enumerate(ref_keys)
        ax     = p_grid[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(θ_circ), S_max .* sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=0.8, label=false)   # thinner circle

        for s in eachindex(oe_lists)
            oe = maps[s][key]
            if length(oe.P_boundary) >= 3
                cx  = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
                idx = sortperm(atan.(oe.Q_boundary .- cy, oe.P_boundary .- cx))
                Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
                Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale

                plot!(ax, Ps, Qs;
                      seriestype = :shape,
                      fillcolor  = colors[s], fillalpha = fillalphas[s],
                      lc = colors[s], lw = 1.2, ls = linestyles[s],  # thinner lw
                      label = false)
            end
        end

        # scatter!(ax, [oe_ref.P_baseline*scale], [oe_ref.Q_baseline*scale];
        #         marker=:star5, ms=5, mc=:crimson, msc=:crimson, label=false)  # smaller star

        hline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)
        vline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)

        # ── Compact title: drop "Node"/"Ph" text to save horizontal space ───
        title!(ax, "N$(oe_ref.unit.node)·$(uppercase(string(oe_ref.unit.phase)))")

        # ── Axis labels: edge cells only, but legible fonts ─────────────────
        row = div(k - 1, ncols) + 1
        col = mod(k - 1, ncols) + 1

        xlabel!(ax, row == nrows ? "P (kW)"   : "")
        ylabel!(ax, col == 1     ? "Q (kvar)" : "")

        if col != 1
            plot!(ax, yformatter = _ -> "", ytickfontsize = 0)
        end
        if row != nrows
            plot!(ax, xformatter = _ -> "", xtickfontsize = 0)
        end

        # Asymmetric margins per cell position
        _lm = col == 1     ? 11Plots.mm  : 1Plots.mm
        _rm = 4Plots.mm    # ← increase for ALL cells, not just last col
        _bm = row == nrows ? 9Plots.mm  : 2.5Plots.mm   # ← 5mm even for interior rows
        _tm = row == 1     ? 3Plots.mm  : 2.5Plots.mm   # ← 6mm for rows 2,3 to absorb xtick bleed

        plot!(ax,
              legend         = false,
              titlefont      = font(12, font_name, :bold),
              tickfont       = font(10, font_name),
              guidefont      = font(11, font_name),
              aspect_ratio   = :equal,
              framestyle     = :box,
              grid           = true,
              gridalpha      = 0.25,
              tick_direction = :in,             # ← :in keeps ticks inside, no external bleed
              ytick_direction = :in,
              yformatter     = col == 1 ? (v -> string(Int(v))) : (_ -> ""),
              xticks         = [-6, -3, 0, 3, 6],
              yticks         = [-6, -3, 0, 3, 6],
              yflip          = false,
              ymirror        = false,           # ← stops right-side y-tick duplication
              xmirror        = false,
              left_margin    = _lm,
              right_margin   = _rm,
              top_margin     = _tm,
              bottom_margin  = _bm)
    end

    for k in n+1 : nrows*ncols
        plot!(p_grid[k], axis=false, ticks=false, grid=false, label=false)
    end

    # ── legend strip (unchanged logic, just full width) ──────────────────────
    p_leg = plot(size             = (W, 55),
                 margin           = 0Plots.mm,
                 background_color = :white,
                 dpi              = dpi,
                 fontfamily       = font_name)

    for s in eachindex(oe_lists)
        plot!(p_leg, Float64[], Float64[];
              seriestype = :shape,
              fillcolor  = colors[s], fillalpha = fillalphas[s],
              lc = colors[s], lw = 1.8, ls = linestyles[s],
              label = names[s])
    end

    plot!(p_leg,
          axis                    = false,
          ticks                   = false,
          grid                    = false,
          legend                  = :inside,
          legendfont              = font(10, font_name),
          legend_columns          = length(oe_lists),
          background_color_inside = :white,
          foreground_color_legend = :black)

    return p_grid, p_leg
end

# helper
fill_img(c, h, w) = fill(c, h, w)


p_grid, p_leg = plot_operating_envelopes_overlay_legend(
                [oe_list11, oe_list12, oe_list21, oe_list22];
                names      = ["M1 baseline", "M1 +VUF", "M2 baseline", "M2 +VUF"],
                colors     = [colorant"#56B4E9", colorant"#0072B2",
                colorant"#E69F00", colorant"#D55E00"],
                linestyles = [:solid, :dash, :solid, :dash],
                fillalphas = [0.2, 0.05, 0.2, 0.05]
)

# display in REPL / VSCode
display(p_grid)
display(p_leg)

# save — two options:

# Option A: save as PNG (raster, fine for submission if dpi=900)
savefig(p_grid, "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_grid.svg")
savefig(p_leg,  "Graphics/overlay_M1_M2_states_TEST$(TEST_NUM)_legend.svg")

# Option B: true single stacked PDF via pdfunite (if you have poppler installed)
# run(`pdfunite Graphics/overlay_M1_M2_states_grid.pdf Graphics/overlay_M1_M2_states_legend.pdf Graphics/overlay_M1_M2_states.pdf`)
