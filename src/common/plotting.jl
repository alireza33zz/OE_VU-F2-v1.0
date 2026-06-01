# ============================================================
# PLOTTING AND VISUALIZATION UTILITIES
# ============================================================
# Color schemes, theme constants, and comprehensive plotting functions
# for three-phase power system visualization.

using Plots, Printf, Colors
gr()

include("structs.jl")
include("network_utils.jl")

# ─────────────────────────────────────────────────────────────────────────────
# THEME CONSTANTS — WHITE BACKGROUND STYLE
# ─────────────────────────────────────────────────────────────────────────────

# Canvas and backgrounds
const _BG       = RGB(1, 1, 1)               # figure canvas (white)
const _BG2      = RGB(1, 1, 1)               # panel surface (white)
const _GRID     = RGBA(0,0,0,0.12)           # light grey grid

# Text and foreground
const _FG       = RGB(0.10, 0.12, 0.16)      # primary text (dark grey)
const _FG2      = RGB(0.35, 0.37, 0.43)      # secondary axis labels (medium grey)

# Phase colors (saturated, high contrast on white)
const _PA  = RGB(0.937, 0.380, 0.216)        # #EF6137  coral-red   (Phase A)
const _PB  = RGB(0.149, 0.737, 0.557)        # #26BC8E  mint-teal   (Phase B)
const _PC  = RGB(0.290, 0.604, 0.976)        # #4A9AF9  sky-blue    (Phase C)

# Pale versions for fills (semi-transparent)
const _PA2 = RGBA(0.937, 0.380, 0.216, 0.28)
const _PB2 = RGBA(0.149, 0.737, 0.557, 0.28)
const _PC2 = RGBA(0.290, 0.604, 0.976, 0.28)

# Demand colors (dark, saturated for emphasis)
const _PA_D = RGB(0.87, 0.30, 0.09)   # deep orange (Phase A demand)
const _PB_D = RGB(0.05, 0.60, 0.45)   # teal-green  (Phase B demand)
const _PC_D = RGB(0.16, 0.47, 0.86)   # sky-blue    (Phase C demand)

# PV fills (very light tints, near-white)
const _PA_PV_FILL = RGB(0.99, 0.94, 0.90)  # very light orange
const _PB_PV_FILL = RGB(0.90, 0.97, 0.94)  # very light teal
const _PC_PV_FILL = RGB(0.92, 0.96, 1.00)  # very light blue

# Hatch pattern for PV bars (works best on GR backend)
const _PV_HATCH = :x

# Accent colors
const _VLIM = RGBA(0.98, 0.58, 0.18, 0.2)   # voltage limit band (orange, semi-transparent)
const _VUF_HI = RGB(0.98, 0.05, 0.05)       # VUF violation indicator (red)

# Common figure configuration kwargs
function _FIG_KW(W::Int, H::Int)
    return (
        size                    = (W, H),
        background_color        = _BG,
        background_color_inside = _BG2,
        foreground_color        = _FG,
        foreground_color_axis   = _FG2,
        foreground_color_border = _FG2,
        foreground_color_grid   = _GRID,
        foreground_color_text   = _FG,
        foreground_color_guide  = _FG,
        dpi                     = 160,
        grid                    = true,
        gridalpha               = 0.14,
        gridlinewidth           = 0.6,
        minorgrid               = false,
        framestyle              = :box,
        tick_direction          = :out,
        fontfamily              = "DejaVu Sans",
    )
end

# Helper: draw a thin zero-line at y=0 (dark grey on white background)
_zeroline!(p) = hline!(p, [0.0]; lc=RGBA(0,0,0,0.25), lw=0.8, label=false)

"""
    _fig_active_power(sys, res) -> Plots.Plot

Visualize three-phase active power (demand and PV injection) per node.

Layout: Three-column cluster per node (phases A, B, C). Within each phase:
  - Downward bar: demand (solid dark color)
  - Upward bar: PV injection (light tint + hatch)

y-axis: Power [kW] (negative = demand, positive = injection)
"""
function _fig_active_power(sys::SystemData, res::OPFResult)
    nodes = collect(1:sys.nN)
    scale = 1e-3   # Convert W to kW

    # Build power matrices (nN × 3)
    Pd = hcat(sys.P_Demand_a, sys.P_Demand_b, sys.P_Demand_c) .* scale
    Pp = hcat(res.P_sol_a,    res.P_sol_b,    res.P_sol_c)    .* scale

    # Clamp near-zero values to exactly zero for clean rendering
    Pd[abs.(Pd) .< 1e-2] .= 0.0
    Pp[abs.(Pp) .< 1e-2] .= 0.0

    # Color schemes
    labels_dem = ["Demand Ph-A", "Demand Ph-B", "Demand Ph-C"]
    labels_pv  = ["PV Ph-A",     "PV Ph-B",     "PV Ph-C"]

    # Cluster offsets: groups three phases with small gaps
    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    δ_inner   = 0.09    # Gap between demand and PV within a phase
    w_bar     = 0.16    # Bar width

    p = plot(; _FIG_KW(1320, 520)...,
        title        = "Active Power — Demand (↓) and PV Injection (↑)",
        titlefont    = (11, _FG, :bold),
        xlabel       = "Node",
        ylabel       = "P  [kW]",
        xticks       = (nodes, string.(nodes)),
        xrotation    = 0,
        ylims        = (-maximum(abs.(Pd))*1.1, maximum(Pp)*1.1),
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legend_columns = 6,
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    # Draw bars for each phase
    for ph in 1:3
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        # Demand bars (downward, saturated)
        bar!(p, x_base .- δ_inner, -Pd[:, ph];
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            label     = labels_dem[ph])

        # PV bars (upward, light fill + hatch + outline)
        bar!(p, x_base .+ δ_inner,  Pp[:, ph];
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,
            linewidth   = 1.1,
            fillpattern = _PV_HATCH,
            label       = labels_pv[ph])
    end

    return p
end

"""
    _fig_reactive_power(sys, res) -> Plots.Plot

Visualize three-phase reactive power (demand and PV injection) per node.
Same layout as _fig_active_power but for Q (reactive power) in kvar.
"""
function _fig_reactive_power(sys::SystemData, res::OPFResult)
    nodes = collect(1:sys.nN)
    scale = 1e-3

    Qd = hcat(sys.Q_Demand_a, sys.Q_Demand_b, sys.Q_Demand_c) .* scale
    Qp = hcat(res.Q_sol_a,    res.Q_sol_b,    res.Q_sol_c)    .* scale

    Qd[abs.(Qd) .< 1e-4] .= 0.0
    Qp[abs.(Qp) .< 1e-4] .= 0.0

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
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legend_columns = 6,
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    for ph in 1:3
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        bar!(p, x_base .- δ_inner, -Qd[:, ph];
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            label     = labels_dem[ph])

        bar!(p, x_base .+ δ_inner,  Qp[:, ph];
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,
            linewidth   = 1.1,
            fillpattern = _PV_HATCH,
            label       = labels_pv[ph])
    end

    return p
end

"""
    _fig_voltage(sys, res, vuf_threshold) -> Plots.Plot

Visualize voltage magnitudes and unbalance factor per node.

Two vertically stacked subplots:
  - Top    : Phase voltages |Va|, |Vb|, |Vc| [V] with ±3% band
  - Bottom : VUF per node [%] as bars (orange if above threshold)
"""
function _fig_voltage(sys::SystemData, res::OPFResult,
                      vuf_threshold::Float64 = 0.02)

    nodes = collect(1:sys.nN)
    vuf   = compute_vuf(res, sys.nN)

    # Convert rectangular (r,i) to magnitude
    _mag(r_, i_) = sqrt.(r_.^2 .+ i_.^2)
    Va = _mag(res.V_r_a, res.V_i_a)
    Vb = _mag(res.V_r_b, res.V_i_b)
    Vc = _mag(res.V_r_c, res.V_i_c)

    Vmin = sys.V_lb[1];  Vmax = sys.V_ub[1]
    Vnom = 400.0 / sqrt(3)

    # Two subplots with different heights
    p = plot(layout  = grid(2, 1; heights=[0.62, 0.38]),
             ; _FIG_KW(1320, 640)...,
             link    = :x,
             bottom_margin = 4Plots.mm,
             left_margin   = 6Plots.mm)

    # ── TOP: VOLTAGE MAGNITUDES ──────────────────────────────────────────────
    
    # Voltage limit band (filled region)
    plot!(p[1],
        [nodes[1]-0.5; nodes; nodes[end]+0.5;
         nodes[end]+0.5; reverse(nodes); nodes[1]-0.5],
        [fill(Vmax, sys.nN+2); fill(Vmin, sys.nN+2)];
        seriestype = :shape,
        fillcolor  = _VLIM,
        linecolor  = :transparent,
        label      = "Standard band")

    # Voltage limit lines (dashed)
    hline!(p[1], [Vmin]; lc=RGBA(0.98,0.01,0.01,0.9), lw=0.9,
           ls=:dash, label=@sprintf("Vmin = %.1f V", Vmin))
    hline!(p[1], [Vmax]; lc=RGBA(0.98,0.01,0.01,0.9), lw=0.9,
           ls=:dash, label=@sprintf("Vmax = %.1f V", Vmax))
    hline!(p[1], [Vnom]; lc=RGBA(0,0,0,0.75), lw=0.8, ls=:dash, label=false)

    # Phase voltage curves
    for (V_, col, lbl) in zip((Va, Vb, Vc), (_PA, _PB, _PC),
                               ("Phase A", "Phase B", "Phase C"))
        plot!(p[1], nodes, V_;
              lc=col, lw=2.0, marker=:circle, ms=4, mc=col, msc=col,
              label=lbl)
    end

    plot!(p[1];
        title      = "Nodal Voltage Magnitudes and Unbalance Factor",
        titlefont  = (11, _FG, :bold),
        ylabel     = "|V|  [V]",
        xlabel     = "",
        xticks     = (nodes, string.(nodes)),
        legend     = :bottomleft,
        legendfont = (7, _FG2),
        ylims      = (Vmin - 1.5, Vmax + 1.5),
    )

    # ── BOTTOM: VUF BARS ─────────────────────────────────────────────────────
    
    thr_pct = vuf_threshold * 100.0
    # Color bars orange if VUF exceeds threshold, else blue
    bar_cols = [v > thr_pct + 0.001 ? RGBA(0.98, 0.55, 0.18, 0.7) : RGBA(0.40, 0.55, 0.80, 0.70) for v in vuf]

    bar!(p[2], nodes, vuf;
         bar_width = 0.6,
         fillcolor = bar_cols,
         linecolor = bar_cols,
         linewidth = 0.3,
         label     = false)

    # VUF threshold line
    hline!(p[2], [thr_pct]; lc=_VUF_HI, lw=1.5, ls=:dash,
           label=@sprintf("Limit %.0f%%", thr_pct))

    plot!(p[2];
        ylabel     = "VUF  [%]",
        xlabel     = "Node",
        xticks     = (nodes, string.(nodes)),
        legend     = :topright,
        legendfont = (7, _FG2),
        ylims      = (0, max(maximum(vuf)*1.35, thr_pct*2.0)),
    )

    return p
end

"""
    _fig_line_flows(sys, res) -> Plots.Plot

Visualize apparent power flows and losses per line.

Two vertically stacked subplots:
  - Top    : |S_ij| per line per phase [kVA] (3-column cluster)
  - Bottom : Active losses per line [W]
"""
function _fig_line_flows(sys::SystemData, res::OPFResult)
    nL = sys.nL
    Pij, Pji, Qij, Qji = compute_line_flows(sys, res)

    # Apparent power at sending end [kVA]
    S_kva = sqrt.(Pij.^2 .+ Qij.^2) ./ 1000.0

    # Losses: sum of |P_ij| - |P_ji| per phase (always ≥ 0)
    losses_W = [sum(abs.(abs.(Pij[k, :]) .- abs.(Pji[k, :]))) for k in 1:nL]

    # Line labels: "i→j"
    xlabels = [@sprintf("%d→%d", sys.fbus[k], sys.tbus[k]) for k in 1:nL]
    xpos    = collect(1:nL)

    phase_lbls = ["Ph-A", "Ph-B", "Ph-C"]

    p = plot(layout  = grid(2, 1; heights=[0.60, 0.40]),
             ; _FIG_KW(1320, 620)...,
             link    = :x,
             bottom_margin = 8Plots.mm,
             left_margin   = 6Plots.mm)

    # ── TOP: |S| PER PHASE ─────────────────────────────────────────────────────
    
    δ_cluster = 0.22
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    w_bar     = 0.18

    for ph in 1:3
        cols_pale  = [_PA2, _PB2, _PC2]
        cols_solid = [_PA, _PB, _PC]
        
        # Hatch if active power direction is negative
        hatch_vec = [Pij[k, ph] < 0 ? :x : :solid for k in 1:nL]

        x_pos = xpos .+ offs_ph[ph]

        bar!(p[1], x_pos, S_kva[:, ph];
             bar_width    = w_bar,
             fillcolor    = cols_pale[ph],
             linecolor    = cols_solid[ph],
             linewidth    = 0.9,
             fillpattern  = hatch_vec,
             label        = phase_lbls[ph])
    end

    # Total apparent power envelope
    S_tot = [sqrt(sum(Pij[k,:].^2 .+ Qij[k,:].^2)) / 1000.0 for k in 1:nL]
    plot!(p[1], xpos, S_tot;
          lc=RGBA(0,0,0,0.55), lw=1.2, ls=:dot,
          marker=:diamond, ms=3, mc=RGBA(0,0,0,0.65), msc=:transparent,
          label="|S| total")

    plot!(p[1];
        title      = "Apparent Power Flows and Active Losses per Line",
        titlefont  = (11, _FG, :bold),
        ylabel     = "|S_ij|  [kVA]",
        xlabel     = "",
        xticks     = (xpos, xlabels),
        xrotation  = 35,
        legend     = :topright,
        legendfont = (8, _FG2),
        ylims      = (0, max(maximum(S_tot) * 1.30, 0.05)),
    )

    # ── BOTTOM: LOSSES ──────────────────────────────────────────────────────────
    
    bar!(p[2], xpos, losses_W;
         bar_width = 0.55,
         fillcolor = RGBA(0.98,0.55,0.18,0.85),
         linecolor = RGBA(0.98,0.55,0.18,0.85),
         linewidth = 0.3,
         label     = false)

    plot!(p[2];
        ylabel     = "ΔP losses  [W]",
        xlabel     = "Line  (from → to)",
        xticks     = (xpos, xlabels),
        xrotation  = 35,
        legend     = false,
        ylims      = (0, maximum(losses_W) * 1.35 + 1e-6),
    )

    return p
end

"""
    plot_dashboard(sys, res; vuf_threshold, save_dir, fmt) -> Dict{Symbol, Plots.Plot}

Generate comprehensive dashboard with four baseline-OPF figures.

Figures:
  - :P_power   : Active power demand + PV injection
  - :Q_power   : Reactive power demand + PV injection
  - :voltage   : Voltage magnitudes + VUF
  - :lines     : Apparent power flows + losses

Parameters:
- save_dir : Directory to save PNG files (empty string = no save)
- fmt      : File extension (default "png")

Returns a Dict allowing individual figure access by key.
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

"""
    plot_operating_envelopes(oe_list) -> Plots.Plot

Plot operating envelopes for all PV units grouped by phase.

Each unit displays its P-Q boundary as a filled region (phase-colored polygon).
Baseline point marked with a small circle.
"""
function plot_operating_envelopes(oe_list::Vector{OEBoundary})
    colors_by_phase = Dict(:a => _PA, :b => _PB, :c => _PC)
    fills_by_phase  = Dict(:a => _PA_PV_FILL, :b => _PB_PV_FILL, :c => _PC_PV_FILL)

    p = plot(; _FIG_KW(1000, 800)...,
        title      = "Operating Envelopes (P-Q Region)",
        titlefont  = (13, _FG, :bold),
        xlabel     = "P  [W]",
        ylabel     = "Q  [var]",
        legend     = :bottomright,
        legendfont = (9, _FG2),
        grid       = true,
    )

    for oe in oe_list
        ph = oe.unit.phase
        col = colors_by_phase[ph]
        fill_col = fills_by_phase[ph]

        # Plot boundary (closed polygon)
        plot!(p, oe.P_boundary, oe.Q_boundary;
            fillalpha  = 0.3,
            fillcolor  = fill_col,
            linecolor  = col,
            linewidth  = 1.5,
            label      = "Node $(oe.unit.node) Ph-$(uppercase(string(ph)))")

        # Mark baseline point
        scatter!(p, [oe.P_baseline], [oe.Q_baseline];
            mc=col, ms=6, msc=col, label=false)
    end

    return p
end
