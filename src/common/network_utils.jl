# ============================================================
# NETWORK UTILITY FUNCTIONS
# ============================================================
# Low-level functions for network matrix assembly, unbalance analysis,
# and line flow calculations.

using LinearAlgebra

include("structs.jl")

"""
    build_admittance_matrix(sys) -> Matrix{ComplexF64}

Assemble the full (3·nN × 3·nN) nodal admittance matrix Y by summing
branch admittances at each node pair, accounting for per-phase coupling.

The admittance matrix encodes the network topology and impedances.
For each branch k with from-bus i and to-bus j:
  - Y[3i-2:3i, 3i-2:3i] += Y_branch  (self-admittance at i)
  - Y[3j-2:3j, 3j-2:3j] += Y_branch  (self-admittance at j)
  - Y[3i-2:3i, 3j-2:3j] -= Y_branch  (mutual admittance i→j)
  - Y[3j-2:3j, 3i-2:3i] -= Y_branch  (mutual admittance j→i)
"""
function build_admittance_matrix(sys::SystemData)::Matrix{ComplexF64}
    n3 = 3 * sys.nN
    Y  = zeros(ComplexF64, n3, n3)

    for k in 1:sys.nL
        Yl = inv(sys.Zij[k])          # Branch admittance (inverse of impedance)
        i = sys.fbus[k];  j = sys.tbus[k]
        ri = (3i-2):(3i);  rj = (3j-2):(3j)
        @inbounds begin
            Y[ri, ri] .+= Yl;  Y[rj, rj] .+= Yl
            Y[ri, rj] .-= Yl;  Y[rj, ri] .-= Yl
        end
    end

    return Y
end

"""
    build_impedance_matrix(Y; slack_node=1) -> Matrix{ComplexF64}

Convert nodal admittance matrix Y to impedance matrix Z by inversion,
excluding the slack bus (typically node 1) which is held at a fixed reference.

Steps:
  1. Remove slack bus rows/columns from Y → Y_red (guaranteed non-singular)
  2. Invert Y_red → Z_red
  3. Embed Z_red back into full-size Z matrix with slack rows/cols = 0

This impedance matrix relates voltage deviations to current injections
at non-slack buses, and is used for sensitivity analysis (e.g., VUF studies).
"""
function build_impedance_matrix(
    Y::Matrix{ComplexF64};
    slack_node::Int = 1
)::Matrix{ComplexF64}

    n3  = size(Y, 1)

    # Indices of the 3 slack bus phases (Node 1 → phases 1,2,3)
    slack_cols = [3*(slack_node-1)+1, 3*(slack_node-1)+2, 3*(slack_node-1)+3]
    keep       = setdiff(1:n3, slack_cols)

    # Reduced admittance matrix (slack bus removed, guaranteed non-singular)
    Y_red = Y[keep, keep]
    Z_red = inv(Y_red)

    # Embed back into full matrix (slack rows/cols remain zero)
    Z_full = zeros(ComplexF64, n3, n3)
    Z_full[keep, keep] = Z_red

    return Z_full
end

"""
    compute_vuf(res, nN) -> Vector{Float64}

Calculate Voltage Unbalance Factor (VUF) at each node [%].

VUF quantifies three-phase voltage asymmetry. Per IEC 61000-4-30:
  VUF = |V₋| / |V₊| × 100%

where V₊ is the positive-sequence (balanced) component and V₋ is negative-sequence.

Algorithm:
  1. For each node, extract phase voltages Va, Vb, Vc (complex)
  2. Decompose into positive and negative sequence using Clarke-Concordia transforms
  3. Compute magnitude of each component
  4. Return ratio as percentage

Returns a vector of VUF values [0, 100+] per node. Typical limits are 1–3%.
"""
function compute_vuf(res::OPFResult, nN::Int)::Vector{Float64}
    # Clarke-Concordia matrix factors: a = e^(j·2π/3)
    ar, ai   = -0.5,  sqrt(3)/2
    a2r, a2i = -0.5, -sqrt(3)/2

    vuf = zeros(nN)
    for i in 1:nN
        Va_r = res.V_r_a[i];  Va_i = res.V_i_a[i]
        Vb_r = res.V_r_b[i];  Vb_i = res.V_i_b[i]
        Vc_r = res.V_r_c[i];  Vc_i = res.V_i_c[i]

        # Positive sequence: Vp = (1/3) · (Va + a·Vb + a²·Vc)
        Vp_r = (1/3)*(Va_r + ar*Vb_r - ai*Vb_i + a2r*Vc_r - a2i*Vc_i)
        Vp_i = (1/3)*(Va_i + ar*Vb_i + ai*Vb_r + a2r*Vc_i + a2i*Vc_r)

        # Negative sequence: Vn = (1/3) · (Va + a²·Vb + a·Vc)
        Vn_r = (1/3)*(Va_r + a2r*Vb_r - a2i*Vb_i + ar*Vc_r - ai*Vc_i)
        Vn_i = (1/3)*(Va_i + a2r*Vb_i + a2i*Vb_r + ar*Vc_i + ai*Vc_r)

        Vp_mag = sqrt(Vp_r^2 + Vp_i^2)
        Vn_mag = sqrt(Vn_r^2 + Vn_i^2)

        # Avoid division by zero for weak grids
        vuf[i] = Vp_mag > 1e-9 ? 100.0 * Vn_mag / Vp_mag : 0.0
    end
    return vuf
end

"""
    compute_line_flows(sys, res) -> (Pij, Pji, Qij, Qji)

Calculate active and reactive power flows on each line per phase.

Returns four nL × 3 matrices:
  - Pij : Active power from bus i to j [W]
  - Pji : Active power from bus j to i [W]
  - Qij : Reactive power from bus i to j [var]
  - Qji : Reactive power from bus j to i [var]

Uses the pi-model with series admittance and half shunt on each end.
Power calculated as: S = V · conj(Y·V) where Y is 6×6 branch admittance.
"""
function compute_line_flows(sys::SystemData, res::OPFResult)
    nL  = sys.nL
    
    # Pre-allocate result matrices
    Pij = zeros(nL, 3);  Pji = zeros(nL, 3)
    Qij = zeros(nL, 3);  Qji = zeros(nL, 3)

    for k in 1:nL
        # Branch admittance matrix (inverse of impedance)
        Yl = inv(sys.Zij[k])
        
        # 6×6 branch admittance matrix (pi-model: series + shunt)
        Yl6 = [(Yl + 0.5 .* sys.Ysh) (-Yl); 
               (-Yl)                (Yl + 0.5 .* sys.Ysh)]

        # Extract from/to bus voltages
        i = sys.fbus[k]
        j = sys.tbus[k]
        
        # Build 6-element complex voltage vector [Vi_a; Vi_b; Vi_c; Vj_a; Vj_b; Vj_c]
        v = [res.V_r_a[i] + im*res.V_i_a[i]; 
             res.V_r_b[i] + im*res.V_i_b[i];
             res.V_r_c[i] + im*res.V_i_c[i]; 
             res.V_r_a[j] + im*res.V_i_a[j];
             res.V_r_b[j] + im*res.V_i_b[j]; 
             res.V_r_c[j] + im*res.V_i_c[j]]

        # Power flow at both ends: S = V · conj(Y·V)
        Sij = v[1:3] .* conj.(Yl6[1:3, :] * v)
        Sji = v[4:6] .* conj.(Yl6[4:6, :] * v)

        # Store real and imaginary parts
        Pij[k, :] = real.(Sij)
        Qij[k, :] = imag.(Sij)
        Pji[k, :] = real.(Sji)
        Qji[k, :] = imag.(Sji)
    end

    return Pij, Pji, Qij, Qji
end
