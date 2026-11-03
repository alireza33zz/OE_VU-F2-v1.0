# ============================================================
# OE PLOTS
# ============================================================

const _θ_CIRC = range(0, 2π, length=200)

# Default 4-scenario overlay style (M1/M2 × baseline/+VUF)
const OVERLAY_NAMES      = ["M1 baseline", "M1 +VUF", "M2 baseline", "M2 +VUF"]
const OVERLAY_COLORS     = [colorant"#56B4E9", colorant"#0072B2", colorant"#E69F00", colorant"#D55E00"]
const OVERLAY_LINESTYLES = [:solid, :dash, :solid, :dash]
const OVERLAY_FILLALPHAS = [0.2, 0.05, 0.2, 0.05]

"""
    _oe_polygon(oe, scale) -> (Ps, Qs)

Closed polygon of the OE boundary, points sorted by angle around their centroid.
"""
function _oe_polygon(oe::OEBoundary, scale)
    cx  = mean(oe.P_boundary);  cy = mean(oe.Q_boundary)
    idx = sortperm(atan.(oe.Q_boundary .- cy, oe.P_boundary .- cx))
    Ps  = [oe.P_boundary[idx]; oe.P_boundary[idx[1]]] .* scale
    Qs  = [oe.Q_boundary[idx]; oe.Q_boundary[idx[1]]] .* scale
    return Ps, Qs
end

"""
    _align_scenarios(oe_lists) -> (ref_keys, maps)

Key every scenario by (node, phase) and check they all cover the reference units.
"""
function _align_scenarios(oe_lists::Vector{Vector{OEBoundary}})
    key_of(oe) = (oe.unit.node, oe.unit.phase)
    ref_keys   = [key_of(oe) for oe in oe_lists[1]]
    maps       = [Dict(key_of(oe) => oe for oe in oelist) for oelist in oe_lists]

    for (i, mp) in enumerate(maps)
        missing_keys = setdiff(ref_keys, collect(Base.keys(mp)))
        isempty(missing_keys) || error("Scenario $i is missing node/phase keys: $missing_keys")
    end
    return ref_keys, maps
end

"""
    plot_operating_envelopes(oe_list; nrows) -> Plots.Plot
"""
function plot_operating_envelopes(oe_list::Vector{OEBoundary};
                                  nrows::Int = floor(Int, sqrt(length(oe_list))))
    n     = length(oe_list)
    ncols = ceil(Int, n / nrows)
    scale = 1e-3   # W → kW

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=900)

    for (k, oe) in enumerate(oe_list)
        ax    = p[k]
        S_max = oe.unit.S_max * scale

        # Hardware S_max circle
        plot!(ax, S_max.*cos.(_θ_CIRC), S_max.*sin.(_θ_CIRC);
              lc=:grey60, ls=:dash, lw=1.2, label="S_max")

        # OE shaded polygon
        if length(oe.P_boundary) >= 3
            Ps, Qs = _oe_polygon(oe, scale)
            plot!(ax, Ps, Qs; seriestype=:shape,
                  fillcolor=:dodgerblue, fillalpha=0.18,
                  lc=:dodgerblue, lw=1.8, label="OE")
        end

        # Baseline operating point
        scatter!(ax, [oe.P_baseline*scale], [oe.Q_baseline*scale];
                 marker=:star5, ms=8, mc=:crimson, msc=:crimson, label="Baseline")

        hline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)
        vline!(ax,[0.0]; lc=:black, lw=0.6, ls=:dot, label=false)

        title!(ax, "Node $(oe.unit.node) · Ph $(uppercase(string(oe.unit.phase)))")
        xlabel!(ax, "P (kW)");  ylabel!(ax, "Q (kvar)")
        plot!(ax, legend=(k==1), legendfont=font(7), titlefont=font(9,:bold),
              tickfont=font(7), guidefont=font(8),
              aspect_ratio=:equal, framestyle=:box, grid=true, gridalpha=0.25)
    end

    # Blank unused grid cells
    for k in n+1:nrows*ncols
        plot!(p[k], axis=false, ticks=false, grid=false, label=false)
    end
    return p
end

"""
    plot_operating_envelopes_overlay(oe_lists; names, colors, fillalphas, linestyles) -> Plot

Overlays multiple OE scenarios (each a Vector{OEBoundary}) on one grid,
aligned by (node, phase). Legend sits in the first panel.
"""
function plot_operating_envelopes_overlay(
    oe_lists::Vector{Vector{OEBoundary}};
    names      = OVERLAY_NAMES,
    colors     = OVERLAY_COLORS,
    fillalphas = OVERLAY_FILLALPHAS,
    linestyles = OVERLAY_LINESTYLES,
    nrows      = 3,
    scale      = 1e-3,
    dpi        = 900
)
    @assert length(oe_lists) == length(names) == length(colors) ==
            length(fillalphas) == length(linestyles)

    ref_keys, maps = _align_scenarios(oe_lists)
    n     = length(ref_keys)
    ncols = ceil(Int, n / nrows)

    p = plot(layout=(nrows,ncols), size=(ncols*320, nrows*320),
             margin=5Plots.mm, background_color=:white, dpi=dpi)

    for (k, key) in enumerate(ref_keys)
        ax     = p[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(_θ_CIRC), S_max .* sin.(_θ_CIRC);
              lc=:grey60, ls=:dash, lw=1.2, label=false)

        for s in eachindex(oe_lists)
            oe = maps[s][key]
            length(oe.P_boundary) >= 3 || continue
            Ps, Qs = _oe_polygon(oe, scale)
            plot!(ax, Ps, Qs;
                  seriestype=:shape,
                  fillcolor=colors[s], fillalpha=fillalphas[s],
                  lc=colors[s], lw=1.8, ls=linestyles[s],
                  label = (k==1 ? names[s] : false))
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

"""
    plot_operating_envelopes_overlay_legend(oe_lists; ..., compact=false) -> (p_grid, p_leg)

Same overlay as `plot_operating_envelopes_overlay`, but without baseline stars and
with the legend returned as a separate full-width strip (for papers).

`compact = true` gives the dense layout used for large networks (IEEE 2-column):
smaller cells, short "N12·A" titles, axis labels/tick labels only on edge cells,
and fixed `ticks` on both axes.
"""
function plot_operating_envelopes_overlay_legend(
    oe_lists::Vector{Vector{OEBoundary}};
    names      = OVERLAY_NAMES,
    colors     = OVERLAY_COLORS,
    fillalphas = OVERLAY_FILLALPHAS,
    linestyles = OVERLAY_LINESTYLES,
    nrows      = 3,
    scale      = 1e-3,
    dpi        = 900,
    font_name  = "Courier",
    compact    = false,
    ticks      = [-6, -3, 0, 3, 6]     # used only when compact = true
)
    @assert length(oe_lists) == length(names) == length(colors) ==
            length(fillalphas) == length(linestyles)

    ref_keys, maps = _align_scenarios(oe_lists)
    n     = length(ref_keys)
    ncols = ceil(Int, n / nrows)
    cell  = compact ? 260 : 320
    W     = ncols * cell
    H     = nrows * cell

    # Line widths: (circle, envelope, zero-lines)
    lw_circ, lw_oe, lw_zero = compact ? (0.8, 1.2, 0.5) : (1.2, 1.8, 0.6)

    # ── main grid (no legend) ────────────────────────────────────────────────
    p_grid = plot(layout = (nrows, ncols),
                  size   = (W, H),
                  margin = compact ? 0Plots.mm : 5Plots.mm,   # compact: margins set per cell
                  background_color = :white,
                  dpi = dpi,
                  fontfamily = font_name)

    for (k, key) in enumerate(ref_keys)
        ax     = p_grid[k]
        oe_ref = maps[1][key]
        S_max  = oe_ref.unit.S_max * scale

        plot!(ax, S_max .* cos.(_θ_CIRC), S_max .* sin.(_θ_CIRC);
              lc=:grey60, ls=:dash, lw=lw_circ, label=false)

        for s in eachindex(oe_lists)
            oe = maps[s][key]
            length(oe.P_boundary) >= 3 || continue
            Ps, Qs = _oe_polygon(oe, scale)
            plot!(ax, Ps, Qs;
                  seriestype = :shape,
                  fillcolor  = colors[s], fillalpha = fillalphas[s],
                  lc = colors[s], lw = lw_oe, ls = linestyles[s],
                  label = false)
        end

        hline!(ax, [0.0]; lc=:black, lw=lw_zero, ls=:dot, label=false)
        vline!(ax, [0.0]; lc=:black, lw=lw_zero, ls=:dot, label=false)

        node_str  = oe_ref.unit.node
        phase_str = uppercase(string(oe_ref.unit.phase))

        if !compact
            title!(ax, "Node $node_str · Ph $phase_str")
            xlabel!(ax, "P (kW)"); ylabel!(ax, "Q (kvar)")
            plot!(ax,
                  legend       = false,
                  titlefont    = font(15, font_name, :bold),
                  tickfont     = font(8, font_name),
                  guidefont    = font(12, font_name),
                  aspect_ratio = :equal,
                  framestyle   = :box,
                  grid         = true,
                  gridalpha    = 0.25)
            continue
        end

        # ── compact layout ───────────────────────────────────────────────────
        title!(ax, "N$(node_str)·$(phase_str)")

        # Axis labels and tick labels on edge cells only
        row = div(k - 1, ncols) + 1
        col = mod(k - 1, ncols) + 1

        xlabel!(ax, row == nrows ? "P (kW)"   : "")
        ylabel!(ax, col == 1     ? "Q (kvar)" : "")

        col != 1     && plot!(ax, yformatter = _ -> "", ytickfontsize = 0)
        row != nrows && plot!(ax, xformatter = _ -> "", xtickfontsize = 0)

        plot!(ax,
              legend          = false,
              titlefont       = font(12, font_name, :bold),
              tickfont        = font(10, font_name),
              guidefont       = font(11, font_name),
              aspect_ratio    = :equal,
              framestyle      = :box,
              grid            = true,
              gridalpha       = 0.25,
              tick_direction  = :in,            # keeps ticks inside, no external bleed
              ytick_direction = :in,
              yformatter      = col == 1 ? (v -> string(Int(v))) : (_ -> ""),
              xticks          = ticks,
              yticks          = ticks,
              yflip           = false,
              ymirror         = false,          # stops right-side y-tick duplication
              xmirror         = false,
              left_margin     = col == 1     ? 11Plots.mm : 1Plots.mm,
              right_margin    = 4Plots.mm,
              top_margin      = row == 1     ? 3Plots.mm  : 2.5Plots.mm,
              bottom_margin   = row == nrows ? 9Plots.mm  : 2.5Plots.mm)
    end

    for k in n+1 : nrows*ncols
        plot!(p_grid[k], axis=false, ticks=false, grid=false, label=false)
    end

    # ── legend-only strip (full width) ───────────────────────────────────────
    p_leg = plot(size             = (W, compact ? 55 : 60),
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
          legendfont              = font(compact ? 10 : 12, font_name),
          legend_columns          = length(oe_lists),
          background_color_inside = :white,
          foreground_color_legend = :black)

    return p_grid, p_leg
end
