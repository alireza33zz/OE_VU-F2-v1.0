# ============================================================
# ADMITTANCE / IMPEDANCE MATRICES
# ============================================================

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

Z = Y⁻¹ with the slack bus removed (slack rows/cols are left at zero).
"""
function build_impedance_matrix(
    Y::Matrix{ComplexF64};
    slack_node::Int = 1
)::Matrix{ComplexF64}

    n3 = size(Y, 1)

    # Indices of the 3 slack bus phases to remove
    slack_cols = [3*(slack_node-1)+1, 3*(slack_node-1)+2, 3*(slack_node-1)+3]
    keep       = setdiff(1:n3, slack_cols)

    # Reduced admittance matrix — guaranteed non-singular
    Z_red = inv(Y[keep, keep])

    # Embed back into full matrix (slack rows/cols remain zero)
    Z_full = zeros(ComplexF64, n3, n3)
    Z_full[keep, keep] = Z_red

    return Z_full
end

"""
    export_admittance_xlsx(Y, filename)

Write real and imaginary parts of Y to two sheets of an Excel file.
"""
function export_admittance_xlsx(Y::Matrix{ComplexF64}, filename::String)
    XLSX.openxlsx(filename, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Matrix_Y_Real")
        XLSX.writetable!(sheet, collect(eachcol(real.(Y))), [Symbol("Col_$i") for i in 1:size(Y,2)])

        sheet2 = XLSX.addsheet!(xf, "Matrix_Y_Imag")
        XLSX.writetable!(sheet2, collect(eachcol(imag.(Y))), [Symbol("Col_$i") for i in 1:size(Y,2)])
    end
    println("  Admittance matrix saved → $filename")
end
