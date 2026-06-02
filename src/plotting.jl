# ============================================================
# plotting.jl — All visualisation helpers
# ============================================================

using Printf

# ─────────────────────────────────────────────────────────────────────────────
# SHARED STYLE CONSTANTS  —  WHITE THEME
# ─────────────────────────────────────────────────────────────────────────────

const _BG       = RGB(1, 1, 1)
const _BG2      = RGB(1, 1, 1)
const _GRID     = RGBA(0,0,0,0.12)
const _FG       = RGB(0.10, 0.12, 0.16)
const _FG2      = RGB(0.35, 0.37, 0.43)

const _PA  = RGB(0.937, 0.380, 0.216)
const _PB  = RGB(0.149, 0.737, 0.557)
const _PC  = RGB(0.290, 0.604, 0.976)

const _PA2 = RGBA(0.937, 0.380, 0.216, 0.28)
const _PB2 = RGBA(0.149, 0.737, 0.557, 0.28)
const _PC2 = RGBA(0.290, 0.604, 0.976, 0.28)

const _PA_D = RGB(0.87, 0.30, 0.09)
const _PB_D = RGB(0.05, 0.60, 0.45)
const _PC_D = RGB(0.16, 0.47, 0.86)

const _PA_PV_FILL = RGB(0.99, 0.94, 0.90)
const _PB_PV_FILL = RGB(0.90, 0.97, 0.94)
const _PC_PV_FILL = RGB(0.92, 0.96, 1.00)

const _PV_HATCH = :x
const _VLIM     = RGBA(0.98, 0.58, 0.18, 0.2)
const _VUF_HI   = RGB(0.98, 0.05, 0.05)

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

_zeroline!(p) = hline!(p, [0.0]; lc=RGBA(0,0,0,0.25), lw=0.8, label=false)

const _NODES = 1:13

const font_name = "Courier"

# ─────────────────────────────────────────────────────────────────────────────
# RESULT REPORTING
# ─────────────────────────────────────────────────────────────────────────────

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

    nN = sys.nN

    mag(r,i) = sqrt.(r.^2 .+ i.^2)
    arg(r,i) = rad2deg.(atan.(i, r))

    Va_mag = mag(res.V_r_a, res.V_i_a);  Va_arg = arg(res.V_r_a, res.V_i_a)
    Vb_mag = mag(res.V_r_b, res.V_i_b);  Vb_arg = arg(res.V_r_b, res.V_i_b)
    Vc_mag = mag(res.V_r_c, res.V_i_c);  Vc_arg = arg(res.V_r_c, res.V_i_c)

    vuf  = compute_vuf(res, nN)
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
    println(DataFrame(
        Line  = 1:sys.nL, From = sys.fbus, To = sys.tbus,
        Pij_a = round.(Pij[:,1]; digits=2), Pij_b = round.(Pij[:,2]; digits=2),
        Pij_c = round.(Pij[:,3]; digits=2),
        Pji_a = round.(Pji[:,1]; digits=2), Pji_b = round.(Pji[:,2]; digits=2),
        Pji_c = round.(Pji[:,3]; digits=2)))

    println("\n── Reactive power flows (var) ──")
    println(DataFrame(
        Line  = 1:sys.nL, From = sys.fbus, To = sys.tbus,
        Qij_a = round.(Qij[:,1]; digits=2), Qij_b = round.(Qij[:,2]; digits=2),
        Qij_c = round.(Qij[:,3]; digits=2),
        Qji_a = round.(Qji[:,1]; digits=2), Qji_b = round.(Qji[:,2]; digits=2),
        Qji_c = round.(Qji[:,3]; digits=2)))
end

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — ACTIVE POWER
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_active_power(sys, res) -> Plots.Plot

Three-column cluster per node (Ph-A, Ph-B, Ph-C). Within each phase, two bars:
  ▼  Active demand  (negative → downward)
  ▲  PV injection   (positive → upward)
"""
function _fig_active_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3

    Pd = hcat(sys.P_Demand_a, sys.P_Demand_b, sys.P_Demand_c) .* scale
    Pp = hcat(res.P_sol_a,    res.P_sol_b,    res.P_sol_c)    .* scale

    Pd[abs.(Pd) .< 1e-2] .= 0.0
    Pp[abs.(Pp) .< 1e-2] .= 0.0

    cols_solid = [_PA,  _PB,  _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Demand Ph-A", "Demand Ph-B", "Demand Ph-C"]
    labels_pv  = ["PV Ph-A",     "PV Ph-B",     "PV Ph-C"]

    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    δ_inner   = 0.09
    w_bar     = 0.16

    ylim_lo = -maximum(abs.(Pd)) * 1.1
    ylim_hi =  maximum(Pp)       * 1.1

    p = plot(; _FIG_KW(1320, 520)...,
        title     = "Active Power — Demand (↓) and PV Injection (↑)",
        titlefont = font(11, :bold),
        xlabel    = "Node",
        ylabel    = "P  [kW]",
        xticks    = (nodes, string.(nodes)),
        xrotation = 0,
        ylims     = (ylim_lo, ylim_hi),
        legend    = :topright)

    for (ph_idx, ph_off) in enumerate(offs_ph)
        xs_dem = nodes .+ ph_off .- δ_inner/2
        xs_pv  = nodes .+ ph_off .+ δ_inner/2

        bar!(p, xs_dem, -abs.(Pd[:, ph_idx]);
             bar_width = w_bar,
             fillcolor = cols_solid[ph_idx], fillalpha = 0.85,
             lc = cols_solid[ph_idx], lw = 0.4,
             label = labels_dem[ph_idx])

        bar!(p, xs_pv, Pp[:, ph_idx];
             bar_width = w_bar,
             fillcolor = cols_pale[ph_idx], fillalpha = 0.9,
             lc = cols_solid[ph_idx], lw = 0.4,
             label = labels_pv[ph_idx])
    end

    _zeroline!(p)
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — REACTIVE POWER
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_reactive_power(sys, res) -> Plots.Plot

Same layout as Fig 1, but for reactive power (kvar).
"""
function _fig_reactive_power(sys::SystemData, res::OPFResult)
    nodes = collect(_NODES)
    scale = 1e-3

    Qd = hcat(sys.Q_Demand_a, sys.Q_Demand_b, sys.Q_Demand_c) .* scale
    Qp = hcat(res.Q_sol_a,    res.Q_sol_b,    res.Q_sol_c)    .* scale

    Qd[abs.(Qd) .< 1e-4] .= 0.0
    Qp[abs.(Qp) .< 1e-4] .= 0.0

    cols_solid = [_PA,  _PB,  _PC]
    cols_pale  = [_PA2, _PB2, _PC2]
    labels_dem = ["Q-demand Ph-A", "Q-demand Ph-B", "Q-demand Ph-C"]
    labels_pv  = ["Q-PV Ph-A",     "Q-PV Ph-B",     "Q-PV Ph-C"]

    δ_cluster = 0.28
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    δ_inner   = 0.09
    w_bar     = 0.16

    all_vals = vcat(vec(-abs.(Qd)), vec(Qp))
    ylim_lo  = minimum(all_vals) * 1.1
    ylim_hi  = maximum(all_vals) * 1.1

    p = plot(; _FIG_KW(1320, 520)...,
        title     = "Reactive Power — Demand (↓) and PV Injection (↑)",
        titlefont = font(11, :bold),
        xlabel    = "Node",
        ylabel    = "Q  [kvar]",
        xticks    = (nodes, string.(nodes)),
        ylims     = (ylim_lo, ylim_hi),
        legend    = :topright)

    for (ph_idx, ph_off) in enumerate(offs_ph)
        xs_dem = nodes .+ ph_off .- δ_inner/2
        xs_pv  = nodes .+ ph_off .+ δ_inner/2

        bar!(p, xs_dem, -abs.(Qd[:, ph_idx]);
             bar_width = w_bar,
             fillcolor = cols_solid[ph_idx], fillalpha = 0.85,
             lc = cols_solid[ph_idx], lw = 0.4,
             label = labels_dem[ph_idx])

        bar!(p, xs_pv, Qp[:, ph_idx];
             bar_width = w_bar,
             fillcolor = cols_pale[ph_idx], fillalpha = 0.9,
             lc = cols_solid[ph_idx], lw = 0.4,
             label = labels_pv[ph_idx])
    end

    _zeroline!(p)
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — VOLTAGE MAGNITUDES + VUF
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_voltage(sys, res, vuf_threshold) -> Plots.Plot

Two vertically stacked subplots:
  Top    : Phase voltage magnitudes |Va|, |Vb|, |Vc| [V] with ±3 % limit band.
  Bottom : VUF [%] bar per node, coloured by severity.
"""
function _fig_voltage(sys::SystemData, res::OPFResult,
                      vuf_threshold::Float64 = 0.02)

    nN    = sys.nN
    nodes = collect(1:nN)
    vuf   = compute_vuf(res, nN)

    _mag(r_, i_) = sqrt.(r_.^2 .+ i_.^2)
    Va = _mag(res.V_r_a, res.V_i_a)
    Vb = _mag(res.V_r_b, res.V_i_b)
    Vc = _mag(res.V_r_c, res.V_i_c)

    Vmin = sys.V_lb[1];  Vmax = sys.V_ub[1]
    Vnom = 400.0 / sqrt(3)

    p = plot(layout        = grid(2, 1; heights=[0.62, 0.38]),
             ; _FIG_KW(1320, 640)...,
             link          = :x,
             bottom_margin = 4Plots.mm,
             left_margin   = 6Plots.mm)

    # ── TOP: voltage magnitudes ──────────────────────────────────────────────
    band_x = [nodes[1]-0.5; nodes; nodes[end]+0.5]
    band_lo = fill(Vmin, length(band_x))
    band_hi = fill(Vmax, length(band_x))

    plot!(p[1], band_x, band_hi;
          fillrange = band_lo,
          fillcolor = _VLIM, fillalpha = 1.0,
          lc = :transparent, label = "V limits")

    plot!(p[1], nodes, Va; lc=_PA, lw=1.8, marker=:circle, ms=4, label="|Va|")
    plot!(p[1], nodes, Vb; lc=_PB, lw=1.8, marker=:circle, ms=4, label="|Vb|")
    plot!(p[1], nodes, Vc; lc=_PC, lw=1.8, marker=:circle, ms=4, label="|Vc|")

    hline!(p[1], [Vnom]; lc=:grey50, lw=0.8, ls=:dot, label="Vnom")

    plot!(p[1],
          title     = "Nodal Voltage Magnitudes",
          titlefont = font(10, :bold),
          ylabel    = "Voltage [V]",
          xlabel    = "",
          xticks    = (nodes, string.(nodes)),
          legend    = :outertopright,
          legendfont = font(8))

    # ── BOTTOM: VUF bars ─────────────────────────────────────────────────────
    bar_colors = [vuf[i] > vuf_threshold*100 ? _VUF_HI : _PB for i in 1:nN]

    bar!(p[2], nodes, vuf;
         bar_width = 0.6,
         fillcolor = bar_colors,
         lc = :transparent,
         label = false)

    hline!(p[2], [vuf_threshold*100]; lc=:red, lw=1.0, ls=:dash, label="Limit")

    plot!(p[2],
          title     = "Voltage Unbalance Factor",
          titlefont = font(10, :bold),
          ylabel    = "VUF [%]",
          xlabel    = "Node",
          xticks    = (nodes, string.(nodes)),
          legend    = :topright,
          legendfont = font(8))

    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — LINE FLOWS + LOSSES
# ─────────────────────────────────────────────────────────────────────────────

"""
    _fig_line_flows(sys, res) -> Plots.Plot

Two vertically stacked subplots:
  Top    : Apparent power |S_ij| per line per phase [kVA].
  Bottom : Active power losses per line [W].
"""
function _fig_line_flows(sys::SystemData, res::OPFResult)
    nL = sys.nL
    Pij, Pji, Qij, Qji = compute_line_flows(sys, res)

    scale = 1e-3
    Sij = sqrt.(Pij.^2 .+ Qij.^2) .* scale   # nL × 3  [kVA]

    losses = [sum(abs.(Pij[l,:]) .- abs.(Pji[l,:])) for l in 1:nL]

    line_labels = ["$(sys.fbus[l])→$(sys.tbus[l])" for l in 1:nL]
    line_idx    = 1:nL

    δ_cluster = 0.22
    offs_ph   = (-δ_cluster, 0.0, +δ_cluster)
    w_bar     = 0.18

    p = plot(layout        = grid(2, 1; heights=[0.60, 0.40]),
             ; _FIG_KW(1320, 640)...,
             bottom_margin = 4Plots.mm,
             left_margin   = 6Plots.mm)

    cols = [_PA, _PB, _PC]
    ph_labels = ["Ph-A", "Ph-B", "Ph-C"]

    for (ph_idx, ph_off) in enumerate(offs_ph)
        xs = collect(line_idx) .+ ph_off
        bar!(p[1], xs, Sij[:, ph_idx];
             bar_width = w_bar,
             fillcolor = cols[ph_idx], fillalpha = 0.85,
             lc = cols[ph_idx], lw = 0.4,
             label = ph_labels[ph_idx])
    end

    plot!(p[1],
          title     = "Apparent Power Flow |S_ij| [kVA]",
          titlefont = font(10, :bold),
          ylabel    = "S [kVA]",
          xlabel    = "",
          xticks    = (collect(line_idx), line_labels),
          xrotation = 45,
          legend    = :topright,
          legendfont = font(8))

    bar!(p[2], collect(line_idx), losses;
         bar_width = 0.55,
         fillcolor = _FG2, fillalpha = 0.75,
         lc = :transparent, label = false)

    plot!(p[2],
          title     = "Active Power Losses [W]",
          titlefont = font(10, :bold),
          ylabel    = "Loss [W]",
          xlabel    = "Line",
          xticks    = (collect(line_idx), line_labels),
          xrotation = 45,
          legend    = false)

    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# DASHBOARD DRIVER
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_dashboard(sys, res; vuf_threshold, save_dir, fmt) -> Dict{Symbol,Plots.Plot}

Build and return all four baseline-OPF figures.
Keys: :P_power, :Q_power, :voltage, :lines
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

# ─────────────────────────────────────────────────────────────────────────────
# OE OVERLAY — basic (no external legend strip)
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_operating_envelopes_overlay(oe_lists; ...) -> Plots.Plot

Overlays multiple OE scenarios on the same grid layout.
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
    ref_keys   = [key_of(oe) for oe in oe_lists[1]]
    n          = length(ref_keys)
    maps       = [Dict(key_of(oe) => oe for oe in oelist) for oelist in oe_lists]

    for (i, mp) in enumerate(maps)
        missing_keys = setdiff(ref_keys, collect(Base.keys(mp)))
        isempty(missing_keys) || error("Scenario $i missing node/phase keys: $missing_keys")
    end

    nrows = 3
    ncols = ceil(Int, n / nrows)

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=dpi)

    θ_circ = range(0, 2π, length=200)

    for (k, key) in enumerate(ref_keys)
        ax     = p[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(θ_circ), S_max .* sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=1.2, label=false)

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
                      lc = colors[s], lw = 1.8, ls = linestyles[s],
                      label = (k == 1 ? names[s] : false))
            end
        end

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

# ─────────────────────────────────────────────────────────────────────────────
# OE OVERLAY — IEEE formatted with external legend strip
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_operating_envelopes_overlay_legend(oe_lists; ...) -> (Plots.Plot, Plots.Plot)

Returns `(p_grid, p_leg)`: the main grid and a separate narrow legend strip.
Designed for 2-column IEEE figures.
"""
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

    nrows  = 3
    ncols  = ceil(Int, n / nrows)
    cell_w = 260
    cell_h = 260
    W      = ncols * cell_w
    H      = nrows * cell_h

    p_grid = plot(layout          = (nrows, ncols),
                  size            = (W, H),
                  margin          = 0Plots.mm,
                  background_color = :white,
                  dpi             = dpi,
                  fontfamily      = font_name)

    θ_circ = range(0, 2π, length=200)

    for (k, key) in enumerate(ref_keys)
        ax     = p_grid[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(θ_circ), S_max .* sin.(θ_circ);
              lc=:grey60, ls=:dash, lw=0.8, label=false)

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
                      lc = colors[s], lw = 1.2, ls = linestyles[s],
                      label = false)
            end
        end

        hline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)
        vline!(ax, [0.0]; lc=:black, lw=0.5, ls=:dot, label=false)

        title!(ax, "N$(oe_ref.unit.node)·$(uppercase(string(oe_ref.unit.phase)))")

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

        _lm = col == 1     ? 11Plots.mm : 1Plots.mm
        _rm = 4Plots.mm
        _bm = row == nrows ?  9Plots.mm : 2.5Plots.mm
        _tm = row == 1     ?  3Plots.mm : 2.5Plots.mm

        plot!(ax,
              legend          = false,
              titlefont       = font(12, font_name, :bold),
              tickfont        = font(10, font_name),
              guidefont       = font(11, font_name),
              aspect_ratio    = :equal,
              framestyle      = :box,
              grid            = true,
              gridalpha       = 0.25,
              tick_direction  = :in,
              ytick_direction = :in,
              yformatter      = col == 1 ? (v -> string(Int(v))) : (_ -> ""),
              xticks          = [-6, -3, 0, 3, 6],
              yticks          = [-6, -3, 0, 3, 6],
              yflip           = false,
              ymirror         = false,
              xmirror         = false,
              left_margin     = _lm,
              right_margin    = _rm,
              top_margin      = _tm,
              bottom_margin   = _bm)
    end

    for k in n+1:nrows*ncols
        plot!(p_grid[k], axis=false, ticks=false, grid=false, label=false)
    end

    # ── Legend strip ─────────────────────────────────────────────────────────
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