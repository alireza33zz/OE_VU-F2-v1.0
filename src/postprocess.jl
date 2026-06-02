# ============================================================
# postprocess.jl — Post-processing of OPF results
# ============================================================

"""
    compute_vuf(res, nN) -> Vector{Float64}

Compute the Voltage Unbalance Factor (%) at every bus from a solved OPFResult.

VUF(%) = 100 · |V₋| / |V₊|
"""
function compute_vuf(res::OPFResult, nN::Int)::Vector{Float64}
    ar, ai   = -0.5,  sqrt(3)/2
    a2r, a2i = -0.5, -sqrt(3)/2

    vuf = zeros(nN)
    for i in 1:nN
        Va_r = res.V_r_a[i];  Va_i = res.V_i_a[i]
        Vb_r = res.V_r_b[i];  Vb_i = res.V_i_b[i]
        Vc_r = res.V_r_c[i];  Vc_i = res.V_i_c[i]

        # Positive sequence
        Vp_r = (1/3)*(Va_r + ar*Vb_r - ai*Vb_i + a2r*Vc_r - a2i*Vc_i)
        Vp_i = (1/3)*(Va_i + ar*Vb_i + ai*Vb_r + a2r*Vc_i + a2i*Vc_r)

        # Negative sequence
        Vn_r = (1/3)*(Va_r + a2r*Vb_r - a2i*Vb_i + ar*Vc_r - ai*Vc_i)
        Vn_i = (1/3)*(Va_i + a2r*Vb_i + a2i*Vb_r + ar*Vc_i + ai*Vc_r)

        Vp_mag = sqrt(Vp_r^2 + Vp_i^2)
        Vn_mag = sqrt(Vn_r^2 + Vn_i^2)

        vuf[i] = Vp_mag > 1e-9 ? 100.0 * Vn_mag / Vp_mag : 0.0
    end
    return vuf
end

function compute_line_flows(sys::SystemData, res::OPFResult)
    nL  = sys.nL
    
    # Pre-allocate result matrices
    Pij = zeros(nL, 3);  Pji = zeros(nL, 3)
    Qij = zeros(nL, 3);  Qji = zeros(nL, 3)

    for k in 1:nL
        # --- 1. Line-specific Admittance Calculation ---
        # Since Zij is now a Vector, we index it by k before inverting
        Yl = inv(sys.Zij[k])
        
        # Build the 6x6 branch admittance matrix for THIS specific line
        # This handles the pi-model (series admittance + shunt capacitance)
        Yl6 = [(Yl + 0.5 .* sys.Ysh) (-Yl); 
               (-Yl)                (Yl + 0.5 .* sys.Ysh)]

        # --- 2. Voltage Vector Extraction ---
        i = sys.fbus[k]
        j = sys.tbus[k]
        
        # Construct the 6-element complex voltage vector [Vi_abc; Vj_abc]
        v = [res.V_r_a[i] + im*res.V_i_a[i]; 
             res.V_r_b[i] + im*res.V_i_b[i];
             res.V_r_c[i] + im*res.V_i_c[i]; 
             res.V_r_a[j] + im*res.V_i_a[j];
             res.V_r_b[j] + im*res.V_i_b[j]; 
             res.V_r_c[j] + im*res.V_i_c[j]]

        # --- 3. Power Flow Calculation ---
        # S = V * conj(I), where I = Y * V
        Sij = v[1:3] .* conj.(Yl6[1:3, :] * v)
        Sji = v[4:6] .* conj.(Yl6[4:6, :] * v)

        # --- 4. Store Results ---
        Pij[k, :] = real.(Sij)
        Qij[k, :] = imag.(Sij)
        Pji[k, :] = real.(Sji)
        Qji[k, :] = imag.(Sji)
    end

    return Pij, Pji, Qij, Qji
end

# ============================================================
# 9. OPERATING ENVELOPE (model reuse; no unsafe warm starts)
# ============================================================

struct PVUnit
    node::Int
    phase::Symbol
    S_max::Float64
end
