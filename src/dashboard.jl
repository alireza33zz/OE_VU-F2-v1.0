# ============================================================
# BASELINE OPF DASHBOARD  (4 figures)  —  WHITE THEME
# ============================================================
#   figs = plot_dashboard(sys, baseline_res; vuf_threshold=0.02, save_dir="opf_figures")
#   display(figs[:P_power]); display(figs[:Q_power]); display(figs[:voltage]); display(figs[:lines])
# ============================================================

# ─────────────────────────────────────────────────────────────────────────────
# SHARED STYLE CONSTANTS  —  WHITE THEME
# ─────────────────────────────────────────────────────────────────────────────

# White canvas
const _BG       = RGB(1, 1, 1)               # figure canvas
const _BG2      = RGB(1, 1, 1)               # panel surface
const _GRID     = RGBA(0,0,0,0.12)           # light grey grid
const _FG       = RGB(0.10, 0.12, 0.16)      # primary text (dark)
const _FG2      = RGB(0.35, 0.37, 0.43)      # secondary / axis text

# Phase colours (saturated, work on white bg)
const _PA  = RGB(0.937, 0.380, 0.216)        # #EF6137  coral-red   phase A
const _PB  = RGB(0.149, 0.737, 0.557)        # #26BC8E  mint-teal   phase B
const _PC  = RGB(0.290, 0.604, 0.976)        # #4A9AF9  sky-blue    phase C

# Pale versions (used for demand/PV fills)
const _PA2 = RGBA(0.937, 0.380, 0.216, 0.28)
const _PB2 = RGBA(0.149, 0.737, 0.557, 0.28)
const _PC2 = RGBA(0.290, 0.604, 0.976, 0.28)

# Demand (strong, saturated)
const _PA_D = RGB(0.87, 0.30, 0.09)   # deep orange (Phase A)
const _PB_D = RGB(0.05, 0.60, 0.45)   # teal-green  (Phase B)
const _PC_D = RGB(0.16, 0.47, 0.86)   # sky-blue    (Phase C)

# PV fills (near-white tints of the same hues)
const _PA_PV_FILL = RGB(0.99, 0.94, 0.90)  # very light orange
const _PB_PV_FILL = RGB(0.90, 0.97, 0.94)  # very light teal
const _PC_PV_FILL = RGB(0.92, 0.96, 1.00)  # very light blue

# Optional: common hatch for PV (works best on GR backend)
const _PV_HATCH = :x

# Voltage limit band colour (subtle on white)
const _VLIM = RGBA(0.98, 0.58, 0.18, 0.2)

# VUF accent
const _VUF_HI = RGB(0.98, 0.05, 0.05)        # orange — VUF > threshold

# Common figure kwargs
_FIG_KW(W, H) = (
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

# Thin spine helper: draw a horizontal zero-line (dark on white)
_zeroline!(p) = hline!(p, [0.0]; lc=RGBA(0,0,0,0.25), lw=0.8, label=false)



# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — ACTIVE POWER (P demand ↓  /  P injection ↑)
# ─────────────────────────────────────────────────────────────────────────────

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
    nodes = collect(1:sys.nN)
    scale = 1e-3   # W → kW

    # --- build bar matrices (nN × 3) -----------------------------------------
    Pd = hcat(sys.P_Demand_a, sys.P_Demand_b, sys.P_Demand_c) .* scale
    Pp = hcat(res.P_sol_a,    res.P_sol_b,    res.P_sol_c)    .* scale

    # Clamp near-zero to exactly 0 for clean rendering
    Pd[abs.(Pd) .< 1e-2] .= 0.0
    Pp[abs.(Pp) .< 1e-2] .= 0.0

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
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legendfontcolor = _FG2,
        legend_columns = 6,   # 3×2 layout for 6 entries
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    # Demand/PV bars grouped by phase with small inner offset
    for ph in 1:3
        # Pick colors per phase
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        # Demand (solid, saturated)
        bar!(p, x_base .- δ_inner, -Pd[:, ph];   # or -Qd for Fig. 2
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            fillalpha = 1.0,
            label     = labels_dem[ph])

        # PV (light tint + hatch + colored outline)
        bar!(p, x_base .+ δ_inner,  Pp[:, ph];   # or  Qp for Fig. 2
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,    # outline matches the phase’s demand color
            linewidth   = 1.1,
            fillalpha   = 1.0,
            fillpattern = _PV_HATCH,     # makes PV immediately distinct
            label       = labels_pv[ph])
    end

    return p
end


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — REACTIVE POWER (Q demand ↓  /  Q injection ↑)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_reactive_power(sys, res) -> Plots.Plot

Same layout as Fig 1, but for reactive power (kvar).
Three-column cluster per node; legend outside (top center), 3×2 layout.
Solar injection uses a distinct hatch fill.
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
        # Legend requested outside, top center
        legend       = :outertop,
        legendfont   = (8, _FG2),
        legendfontcolor = _FG2,
        legend_columns = 6,
        bottom_margin = 6Plots.mm,
        left_margin   = 6Plots.mm,
    )

    _zeroline!(p)

    for ph in 1:3
        # Phase-based colors (same spectrum; demand strong, PV very light)
        demand_col = ph == 1 ? _PA_D : ph == 2 ? _PB_D : _PC_D
        pv_fill    = ph == 1 ? _PA_PV_FILL : ph == 2 ? _PB_PV_FILL : _PC_PV_FILL

        x_base = nodes .+ offs_ph[ph]

        # Demand (downward, saturated fill)
        bar!(p, x_base .- δ_inner, -Qd[:, ph];
            bar_width = w_bar,
            fillcolor = demand_col,
            linecolor = demand_col,
            linewidth = 0.0,
            fillalpha = 1.0,
            label     = labels_dem[ph])

        # PV (upward, very light tint + colored outline + hatch)
        bar!(p, x_base .+ δ_inner,  Qp[:, ph];
            bar_width   = w_bar,
            fillcolor   = pv_fill,
            linecolor   = demand_col,   # phase-matched outline
            linewidth   = 1.1,
            fillalpha   = 1.0,
            fillpattern = _PV_HATCH,    # e.g., :x (requires GR or a supported backend)
            label       = labels_pv[ph])
    end

    return p
end


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — VOLTAGE MAGNITUDES + VUF
# ─────────────────────────────────────────────────────────────────────────────

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
    nodes = collect(1:sys.nN)
    vuf   = compute_vuf(res, nN)

    _mag(r_, i_) = sqrt.(r_.^2 .+ i_.^2)
    Va = _mag(res.V_r_a, res.V_i_a)
    Vb = _mag(res.V_r_b, res.V_i_b)
    Vc = _mag(res.V_r_c, res.V_i_c)

    Vmin = sys.V_lb[1];  Vmax = sys.V_ub[1]
    Vnom = sys.V_nom

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
         nodes[end]+0.5; reverse(nodes); nodes[1]-0.5],
        [fill(Vmax, nN+2); fill(Vmin, nN+2)];
        seriestype = :shape,
        fillcolor  = _VLIM,
        linecolor  = :transparent,
        label      = "Standard band")

    # Vmin / Vmax dashed lines
    hline!(p[1], [Vmin]; lc=RGBA(0.98,0.01,0.01,0.9), lw=0.9,
           ls=:dash, label=@sprintf("Vmin = %.1f V", Vmin))
    hline!(p[1], [Vmax]; lc=RGBA(0.98,0.01,0.01,0.9), lw=0.9,
           ls=:dash, label=@sprintf("Vmax = %.1f V", Vmax))
    hline!(p[1], [Vnom]; lc=RGBA(0,0,0,0.75), lw=0.8, ls=:dash, label=false)

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
        legendfontcolor = _FG2,
        ylims      = (Vmin - 1.5, Vmax + 1.5),
    )

    # ── BOTTOM: VUF bars ──────────────────────────────────────────────────────
    # Colour each bar: orange if VUF exceeds threshold, else muted blue-grey
    thr_pct = vuf_threshold * 100.0
    bar_cols = [v > thr_pct + 0.001 ? RGBA(0.98, 0.55, 0.18, 0.7) : RGBA(0.40, 0.55, 0.80, 0.70) for v in vuf]

    bar!(p[2], nodes, vuf;
         bar_width = 0.6,
         fillcolor = bar_cols,
         linecolor = bar_cols,
         linewidth = 0.3,
         fillalpha = 1.0,
         label     = false)

    # Threshold line
    hline!(p[2], [thr_pct]; lc=_VUF_HI, lw=1.5, ls=:dash,
           label=@sprintf("Limit %.0f%%", thr_pct))

    plot!(p[2];
        ylabel     = "VUF  [%]",
        xlabel     = "Node",
        xticks     = (nodes, string.(nodes)),
        legend     = :topright,
        legendfont = (7, _FG2),
        legendfontcolor = _FG2,
        ylims      = (0, max(maximum(vuf)*1.35, thr_pct*2.0)),
    )

    return p
end


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — APPARENT POWER FLOWS IN LINES
# ─────────────────────────────────────────────────────────────────────────────

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
    S_kva = sqrt.(Pij.^2 .+ Qij.^2) ./ 1000.0

    # Losses per line — sum over phases of | |P_ij| − |P_ji| |  (always ≥ 0)
    losses_W = [sum(abs.(abs.(Pij[k, :]) .- abs.(Pji[k, :]))) for k in 1:nL]

    # x labels: "i→j"
    xlabels = [@sprintf("%d→%d", sys.fbus[k], sys.tbus[k]) for k in 1:nL]
    xpos    = collect(1:nL)

    cols_solid = [_PA, _PB, _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    phase_lbls = ["Ph-A", "Ph-B", "Ph-C"]

    p = plot(layout  = grid(2, 1; heights=[0.60, 0.40]),
             ; _FIG_KW(1320, 620)...,
             link    = :x,
             bottom_margin = 8Plots.mm,
             left_margin   = 6Plots.mm)

    # ── TOP: |S| per phase per line — 3-column cluster per line ──────────────
    # Cluster offsets for phases around each integer x = line index
    δ_cluster = 0.22
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    w_bar     = 0.18

    for ph in 1:3
        # Hatch any bar where active power direction is negative (P_ij < 0)
        hatch_vec = [Pij[k, ph] < 0 ? :x : :solid for k in 1:nL]   # per-bar

        # Shifted x positions for this phase
        x_pos = xpos .+ offs_ph[ph]

        bar!(p[1], x_pos, S_kva[:, ph];
             bar_width    = w_bar,
             # no :dodge, we physically shift x to form three-column units
             fillcolor    = cols_pale[ph],
             linecolor    = cols_solid[ph],
             linewidth    = 0.9,
             fillalpha    = 1.0,
             fillpattern  = hatch_vec,
             label        = phase_lbls[ph])
    end

    # Total |S| envelope (thin line) at integer positions
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
        legendfontcolor = _FG2,
        ylims      = (0, max(maximum(S_tot) * 1.30, 0.05)),
    )

    # ── BOTTOM: active losses per line (non-negative) ─────────────────────────
    bar!(p[2], xpos, losses_W;
         bar_width = 0.55,
         fillcolor = RGBA(0.98,0.55,0.18,0.85),
         linecolor = RGBA(0.98,0.55,0.18,0.85),
         linewidth = 0.3,
         fillalpha = 1.0,
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


# ─────────────────────────────────────────────────────────────────────────────
# MASTER FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

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
