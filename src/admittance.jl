# ============================================================
# admittance.jl — Nodal admittance and impedance matrices
# ============================================================

"""
    build_admittance_matrix(sys) -> Matrix{ComplexF64}

Assemble the full (3·nN × 3·nN) nodal admittance matrix Y.
Supports a single global impedance matrix or a per-line Vector of matrices.
"""
"""
    build_admittance_matrix(sys) -> Matrix{ComplexF64}

Assemble the full (3·nN × 3·nN) nodal admittance matrix Y by looping over the line list.
"""

function build_admittance_matrix(sys::SystemData)::Matrix{ComplexF64}
    n3 = 3 * sys.nN
    Y  = zeros(ComplexF64, n3, n3)

    for k in 1:sys.nL
        Yl = inv(sys.Zij[k])          # k-th branch matrix
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

Invert the reduced admittance matrix (slack rows/cols removed) and embed
the result back into a full (3·nN × 3·nN) matrix.
"""
function build_impedance_matrix(
    Y::Matrix{ComplexF64};
    slack_node::Int = 1
)::Matrix{ComplexF64}

    n3  = size(Y, 1)

    # Indices of the 3 slack bus phases to remove
    slack_cols = [3*(slack_node-1)+1, 3*(slack_node-1)+2, 3*(slack_node-1)+3]
    keep       = setdiff(1:n3, slack_cols)

    # Reduced admittance matrix — guaranteed non-singular
    Y_red = Y[keep, keep]
    Z_red = inv(Y_red)

    # Embed back into full matrix (slack rows/cols remain zero)
    Z_full = zeros(ComplexF64, n3, n3)
    Z_full[keep, keep] = Z_red

    return Z_full
end
