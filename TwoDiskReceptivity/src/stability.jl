# Stability equation coefficient matrix assembly

using LinearAlgebra

"""
    assemble_coeff_matrix_BEK1(F, G, H, R, N_cheb, D, D2, Res)

Assemble the coefficient matrices for the linear stability equations
using the single-disk Bödewadt scaling (mode 1).

The nondimensionalisation uses √Res as the wall-normal length scale.
"""
function assemble_coeff_matrix_BEK1(F, G, H, R, N_cheb, D, D2, Res)
    Res0 = sqrt(Res)
    size = N_cheb + 1
    eye = I(N_cheb + 1)
    Zero = zeros(N_cheb + 1, N_cheb + 1)

    # Block (1,1): continuity
    Ta_11 = eye
    A_11 = F .* eye
    B_11 = (1/R) * G .* eye
    C_11 = (1/(R*Res0)) * H .* eye
    dC_11 = D * diag(C_11) .* eye  # D acting on C_11
    D_11 = (1/R) * F .* eye
    D_12 = -(1/R) * 2 * (G .+ 1) .* eye
    D_13 = 1/Res0 * D * F .* eye
    Vxx_11 = -(1/R) * eye
    Vyy_11 = -(1/R^3) * eye
    Vzz_11 = -(1/(R*Res)) * eye
    dVzz_11 = D * diag(Vzz_11) .* eye
    d2Vzz_11 = D2 * diag(Vzz_11) .* eye

    # Block (1,2)
    Ta_12 = Zero;    A_12 = Zero;    B_12 = Zero
    C_12 = Zero;     dC_12 = Zero;   D_14 = Zero
    Vxx_12 = Zero;   Vyy_12 = Zero;  Vzz_12 = Zero
    dVzz_12 = Zero;  d2Vzz_12 = Zero

    # Block (1,3)
    Ta_13 = Zero;    A_13 = Zero;    B_13 = Zero
    C_13 = Zero;     dC_13 = Zero
    Vxx_13 = Zero;   Vyy_13 = Zero;  Vzz_13 = Zero
    dVzz_13 = Zero;  d2Vzz_13 = Zero

    # Block (1,4)
    Ta_14 = Zero;    A_14 = eye;     B_14 = Zero
    C_14 = Zero;     dC_14 = Zero
    Vxx_14 = Zero;   Vyy_14 = Zero;  Vzz_14 = Zero
    dVzz_14 = Zero;  d2Vzz_14 = Zero

    # Block (2,1): radial momentum
    Ta_21 = Zero
    A_21 = Zero
    B_21 = Zero
    C_21 = Zero;     dC_21 = Zero
    D_21 = (1/R) * 2 * (G .+ 1) .* eye
    Vxx_21 = Zero;   Vyy_21 = Zero;  Vzz_21 = Zero
    dVzz_21 = Zero;  d2Vzz_21 = Zero

    # Block (2,2)
    Ta_22 = eye
    A_22 = F .* eye
    B_22 = (1/R) * G .* eye
    C_22 = (1/(R*Res0)) * H .* eye
    dC_22 = D * diag(C_22) .* eye
    D_22 = (1/R) * F .* eye
    D_23 = 1/Res0 * D * G .* eye
    Vxx_22 = -(1/R) * eye
    Vyy_22 = -(1/R^3) * eye
    Vzz_22 = -(1/(R*Res)) * eye
    dVzz_22 = D * diag(Vzz_22) .* eye
    d2Vzz_22 = D2 * diag(Vzz_22) .* eye

    # Block (2,3)
    Ta_23 = Zero;    A_23 = Zero;    B_23 = Zero
    C_23 = Zero;     dC_23 = Zero
    Vxx_23 = Zero;   Vyy_23 = Zero;  Vzz_23 = Zero
    dVzz_23 = Zero;  d2Vzz_23 = Zero

    # Block (2,4)
    Ta_24 = Zero;    A_24 = Zero;    B_24 = (1/R) .* eye
    C_24 = Zero;     dC_24 = Zero;   D_24 = Zero
    Vxx_24 = Zero;   Vyy_24 = Zero;  Vzz_24 = Zero
    dVzz_24 = Zero;  d2Vzz_24 = Zero

    # Block (3,1): azimuthal momentum
    Ta_31 = Zero;    A_31 = Zero;    B_31 = Zero
    C_31 = Zero;     dC_31 = Zero;   D_31 = Zero
    Vxx_31 = Zero;   Vyy_31 = Zero;  Vzz_31 = Zero
    dVzz_31 = Zero;  d2Vzz_31 = Zero

    # Block (3,2)
    Ta_32 = Zero;    A_32 = Zero;    B_32 = Zero
    C_32 = Zero;     dC_32 = Zero;   D_32 = Zero
    Vxx_32 = Zero;   Vyy_32 = Zero;  Vzz_32 = Zero
    dVzz_32 = Zero;  d2Vzz_32 = Zero

    # Block (3,3)
    Ta_33 = eye
    A_33 = F .* eye
    B_33 = (1/R) * G .* eye
    C_33 = (1/(R*Res0)) * H .* eye
    dC_33 = D * diag(C_33) .* eye
    C_34 = 1/Res0 * eye
    dC_34 = D * diag(C_34) .* eye
    D_33 = (1/(R*Res0)) * D * H .* eye
    Vxx_33 = -(1/R) * eye
    Vyy_33 = -(1/R^3) * eye
    Vzz_33 = -(1/(R*Res)) * eye
    dVzz_33 = D * diag(Vzz_33) .* eye
    d2Vzz_33 = D2 * diag(Vzz_33) .* eye

    # Block (3,4)
    Ta_34 = Zero;    A_34 = Zero;    B_34 = Zero
    D_34 = Zero
    Vxx_34 = Zero;   Vyy_34 = Zero;  Vzz_34 = Zero
    dVzz_34 = Zero;  d2Vzz_34 = Zero

    # Block (4,1): wall-normal momentum
    Ta_41 = Zero
    A_41 = eye
    B_41 = Zero
    C_41 = Zero;     dC_41 = Zero
    D_41 = 1/R .* eye
    Vxx_41 = Zero;   Vyy_41 = Zero;  Vzz_41 = Zero
    dVzz_41 = Zero;  d2Vzz_41 = Zero

    # Block (4,2)
    Ta_42 = Zero;    A_42 = Zero;    B_42 = (1/R) .* eye
    C_42 = Zero;     dC_42 = Zero;   D_42 = Zero
    Vxx_42 = Zero;   Vyy_42 = Zero;  Vzz_42 = Zero
    dVzz_42 = Zero;  d2Vzz_42 = Zero

    # Block (4,3)
    Ta_43 = Zero;    A_43 = Zero;    B_43 = Zero
    C_43 = (1/Res0) * eye
    dC_43 = D * diag(C_43) .* eye
    D_43 = Zero
    Vxx_43 = Zero;   Vyy_43 = Zero;  Vzz_43 = Zero
    dVzz_43 = Zero;  d2Vzz_43 = Zero

    # Block (4,4)
    Ta_44 = Zero;    A_44 = Zero;    B_44 = Zero
    C_44 = Zero;     dC_44 = Zero;   D_44 = Zero
    Vxx_44 = Zero;   Vyy_44 = Zero;  Vzz_44 = Zero
    dVzz_44 = Zero;  d2Vzz_44 = Zero

    # Assemble full 4×4 block matrices
    Ta = [Ta_11 Ta_12 Ta_13 Ta_14; Ta_21 Ta_22 Ta_23 Ta_24; Ta_31 Ta_32 Ta_33 Ta_34; Ta_41 Ta_42 Ta_43 Ta_44]
    A  = [A_11  A_12  A_13  A_14;  A_21  A_22  A_23  A_24;  A_31  A_32  A_33  A_34;  A_41  A_42  A_43  A_44]
    B  = [B_11  B_12  B_13  B_14;  B_21  B_22  B_23  B_24;  B_31  B_32  B_33  B_34;  B_41  B_42  B_43  B_44]
    Cm = [C_11  C_12  C_13  C_14;  C_21  C_22  C_23  C_24;  C_31  C_32  C_33  C_34;  C_41  C_42  C_43  C_44]
    dC = [dC_11 dC_12 dC_13 dC_14; dC_21 dC_22 dC_23 dC_24; dC_31 dC_32 dC_33 dC_34; dC_41 dC_42 dC_43 dC_44]
    D1 = [D_11  D_12  D_13  D_14;  D_21  D_22  D_23  D_24;  D_31  D_32  D_33  D_34;  D_41  D_42  D_43  D_44]
    Vxx = [Vxx_11 Vxx_12 Vxx_13 Vxx_14; Vxx_21 Vxx_22 Vxx_23 Vxx_24; Vxx_31 Vxx_32 Vxx_33 Vxx_34; Vxx_41 Vxx_42 Vxx_43 Vxx_44]
    Vyy = [Vyy_11 Vyy_12 Vyy_13 Vyy_14; Vyy_21 Vyy_22 Vyy_23 Vyy_24; Vyy_31 Vyy_32 Vyy_33 Vyy_34; Vyy_41 Vyy_42 Vyy_43 Vyy_44]
    Vzz = [Vzz_11 Vzz_12 Vzz_13 Vzz_14; Vzz_21 Vzz_22 Vzz_23 Vzz_24; Vzz_31 Vzz_32 Vzz_33 Vzz_34; Vzz_41 Vzz_42 Vzz_43 Vzz_44]
    dVzz = [dVzz_11 dVzz_12 dVzz_13 dVzz_14; dVzz_21 dVzz_22 dVzz_23 dVzz_24; dVzz_31 dVzz_32 dVzz_33 dVzz_34; dVzz_41 dVzz_42 dVzz_43 dVzz_44]
    d2Vzz = [d2Vzz_11 d2Vzz_12 d2Vzz_13 d2Vzz_14; d2Vzz_21 d2Vzz_22 d2Vzz_23 d2Vzz_24; d2Vzz_31 d2Vzz_32 d2Vzz_33 d2Vzz_34; d2Vzz_41 d2Vzz_42 d2Vzz_43 d2Vzz_44]

    # Unused blocks (Vxy, Vxz, Vyz) – initialised to zero
    Vxy = zeros(4*size, 4*size)
    Vxz = zeros(4*size, 4*size)
    Vyz = zeros(4*size, 4*size)
    dVxz = zeros(4*size, 4*size)
    dVyz = zeros(4*size, 4*size)

    return CoeffMatrix(Ta, A, B, Cm, dC, D1, Vxx, Vyy, Vzz, dVzz, d2Vzz, Vxy, Vxz, dVxz, Vyz, dVyz)
end

"""
    assemble_coeff_matrix_BEK2(F, G, H, R, N_cheb, D, D2, Ro, Co)

Assemble coefficient matrices for the two-disk BEK formulation with rotation ratio Ro.
The scaling differs from BEK1 — uses Ro in the coefficients instead of 1/√Res.
"""
function assemble_coeff_matrix_BEK2(F, G, H, R, N_cheb, D, D2, Ro, Co)
    if Ro ≈ -1.0
        Ro = 1.0  # Handle Ro = -1 case
    end
    size = N_cheb + 1
    eye = I(N_cheb + 1)
    Zero = zeros(N_cheb + 1, N_cheb + 1)

    # Block (1,1): continuity
    Ta_11 = eye
    A_11 = F .* eye
    B_11 = (1/R) * G .* eye
    C_11 = (Ro/R) * H .* eye
    dC_11 = D * diag(C_11) .* eye
    D_11 = (Ro/R) * F .* eye
    D_12 = -(1/R) * (2*Ro*G .+ Co) .* eye
    D_13 = D * F .* eye
    Vxx_11 = -(1/R) * eye
    Vyy_11 = -(1/R^3) * eye
    Vzz_11 = -(1/R) * eye
    dVzz_11 = D * diag(Vzz_11) .* eye
    d2Vzz_11 = D2 * diag(Vzz_11) .* eye

    # Block (1,2)
    Ta_12 = Zero;    A_12 = Zero;    B_12 = Zero
    C_12 = Zero;     dC_12 = Zero;   D_14 = Zero
    Vxx_12 = Zero;   Vyy_12 = Zero;  Vzz_12 = Zero
    dVzz_12 = Zero;  d2Vzz_12 = Zero

    # Block (1,3)
    Ta_13 = Zero;    A_13 = Zero;    B_13 = Zero
    C_13 = Zero;     dC_13 = Zero
    Vxx_13 = Zero;   Vyy_13 = Zero;  Vzz_13 = Zero
    dVzz_13 = Zero;  d2Vzz_13 = Zero

    # Block (1,4)
    Ta_14 = Zero;    A_14 = eye;     B_14 = Zero
    C_14 = Zero;     dC_14 = Zero
    Vxx_14 = Zero;   Vyy_14 = Zero;  Vzz_14 = Zero
    dVzz_14 = Zero;  d2Vzz_14 = Zero

    # Block (2,1): radial momentum
    Ta_21 = Zero
    A_21 = Zero;     B_21 = Zero
    C_21 = Zero;     dC_21 = Zero
    D_21 = (1/R) * (2*Ro*G .+ Co) .* eye
    Vxx_21 = Zero;   Vyy_21 = Zero;  Vzz_21 = Zero
    dVzz_21 = Zero;  d2Vzz_21 = Zero

    # Block (2,2)
    Ta_22 = eye
    A_22 = F .* eye
    B_22 = (1/R) * G .* eye
    C_22 = (Ro/R) * H .* eye
    dC_22 = D * diag(C_22) .* eye
    D_22 = (Ro/R) * F .* eye
    D_23 = D * G .* eye
    Vxx_22 = -(1/R) * eye
    Vyy_22 = -(1/R^3) * eye
    Vzz_22 = -(1/R) * eye
    dVzz_22 = D * diag(Vzz_22) .* eye
    d2Vzz_22 = D2 * diag(Vzz_22) .* eye

    # Block (2,3)
    Ta_23 = Zero;    A_23 = Zero;    B_23 = Zero
    C_23 = Zero;     dC_23 = Zero
    Vxx_23 = Zero;   Vyy_23 = Zero;  Vzz_23 = Zero
    dVzz_23 = Zero;  d2Vzz_23 = Zero

    # Block (2,4)
    Ta_24 = Zero;    A_24 = Zero;    B_24 = (1/R) .* eye
    C_24 = Zero;     dC_24 = Zero;   D_24 = Zero
    Vxx_24 = Zero;   Vyy_24 = Zero;  Vzz_24 = Zero
    dVzz_24 = Zero;  d2Vzz_24 = Zero

    # Block (3,1): azimuthal momentum
    Ta_31 = Zero;    A_31 = Zero;    B_31 = Zero
    C_31 = Zero;     dC_31 = Zero;   D_31 = Zero
    Vxx_31 = Zero;   Vyy_31 = Zero;  Vzz_31 = Zero
    dVzz_31 = Zero;  d2Vzz_31 = Zero

    # Block (3,2)
    Ta_32 = Zero;    A_32 = Zero;    B_32 = Zero
    C_32 = Zero;     dC_32 = Zero;   D_32 = Zero
    Vxx_32 = Zero;   Vyy_32 = Zero;  Vzz_32 = Zero
    dVzz_32 = Zero;  d2Vzz_32 = Zero

    # Block (3,3)
    Ta_33 = eye
    A_33 = F .* eye
    B_33 = (1/R) * G .* eye
    C_33 = (Ro/R) * H .* eye
    dC_33 = D * diag(C_33) .* eye
    C_34 = eye
    dC_34 = D * diag(C_34) .* eye
    D_33 = (Ro/R) * D * H .* eye
    Vxx_33 = -(1/R) * eye
    Vyy_33 = -(1/R^3) * eye
    Vzz_33 = -(1/R) * eye
    dVzz_33 = D * diag(Vzz_33) .* eye
    d2Vzz_33 = D2 * diag(Vzz_33) .* eye

    # Block (3,4)
    Ta_34 = Zero;    A_34 = Zero;    B_34 = Zero
    D_34 = Zero
    Vxx_34 = Zero;   Vyy_34 = Zero;  Vzz_34 = Zero
    dVzz_34 = Zero;  d2Vzz_34 = Zero

    # Block (4,1): wall-normal momentum
    Ta_41 = Zero
    A_41 = eye
    B_41 = Zero
    C_41 = Zero;     dC_41 = Zero
    D_41 = Ro/R .* eye
    Vxx_41 = Zero;   Vyy_41 = Zero;  Vzz_41 = Zero
    dVzz_41 = Zero;  d2Vzz_41 = Zero

    # Block (4,2)
    Ta_42 = Zero;    A_42 = Zero;    B_42 = (1/R) .* eye
    C_42 = Zero;     dC_42 = Zero;   D_42 = Zero
    Vxx_42 = Zero;   Vyy_42 = Zero;  Vzz_42 = Zero
    dVzz_42 = Zero;  d2Vzz_42 = Zero

    # Block (4,3)
    Ta_43 = Zero;    A_43 = Zero;    B_43 = Zero
    C_43 = eye
    dC_43 = D * diag(C_43) .* eye
    D_43 = Zero
    Vxx_43 = Zero;   Vyy_43 = Zero;  Vzz_43 = Zero
    dVzz_43 = Zero;  d2Vzz_43 = Zero

    # Block (4,4)
    Ta_44 = Zero;    A_44 = Zero;    B_44 = Zero
    C_44 = Zero;     dC_44 = Zero;   D_44 = Zero
    Vxx_44 = Zero;   Vyy_44 = Zero;  Vzz_44 = Zero
    dVzz_44 = Zero;  d2Vzz_44 = Zero

    # Assemble full 4×4 block matrices
    Ta = [Ta_11 Ta_12 Ta_13 Ta_14; Ta_21 Ta_22 Ta_23 Ta_24; Ta_31 Ta_32 Ta_33 Ta_34; Ta_41 Ta_42 Ta_43 Ta_44]
    A  = [A_11  A_12  A_13  A_14;  A_21  A_22  A_23  A_24;  A_31  A_32  A_33  A_34;  A_41  A_42  A_43  A_44]
    B  = [B_11  B_12  B_13  B_14;  B_21  B_22  B_23  B_24;  B_31  B_32  B_33  B_34;  B_41  B_42  B_43  B_44]
    Cm = [C_11  C_12  C_13  C_14;  C_21  C_22  C_23  C_24;  C_31  C_32  C_33  C_34;  C_41  C_42  C_43  C_44]
    dC = [dC_11 dC_12 dC_13 dC_14; dC_21 dC_22 dC_23 dC_24; dC_31 dC_32 dC_33 dC_34; dC_41 dC_42 dC_43 dC_44]
    D1 = [D_11  D_12  D_13  D_14;  D_21  D_22  D_23  D_24;  D_31  D_32  D_33  D_34;  D_41  D_42  D_43  D_44]
    Vxx = [Vxx_11 Vxx_12 Vxx_13 Vxx_14; Vxx_21 Vxx_22 Vxx_23 Vxx_24; Vxx_31 Vxx_32 Vxx_33 Vxx_34; Vxx_41 Vxx_42 Vxx_43 Vxx_44]
    Vyy = [Vyy_11 Vyy_12 Vyy_13 Vyy_14; Vyy_21 Vyy_22 Vyy_23 Vyy_24; Vyy_31 Vyy_32 Vyy_33 Vyy_34; Vyy_41 Vyy_42 Vyy_43 Vyy_44]
    Vzz = [Vzz_11 Vzz_12 Vzz_13 Vzz_14; Vzz_21 Vzz_22 Vzz_23 Vzz_24; Vzz_31 Vzz_32 Vzz_33 Vzz_34; Vzz_41 Vzz_42 Vzz_43 Vzz_44]
    dVzz = [dVzz_11 dVzz_12 dVzz_13 dVzz_14; dVzz_21 dVzz_22 dVzz_23 dVzz_24; dVzz_31 dVzz_32 dVzz_33 dVzz_34; dVzz_41 dVzz_42 dVzz_43 dVzz_44]
    d2Vzz = [d2Vzz_11 d2Vzz_12 d2Vzz_13 d2Vzz_14; d2Vzz_21 d2Vzz_22 d2Vzz_23 d2Vzz_24; d2Vzz_31 d2Vzz_32 d2Vzz_33 d2Vzz_34; d2Vzz_41 d2Vzz_42 d2Vzz_43 d2Vzz_44]

    Vxy = zeros(4*size, 4*size)
    Vxz = zeros(4*size, 4*size)
    Vyz = zeros(4*size, 4*size)
    dVxz = zeros(4*size, 4*size)
    dVyz = zeros(4*size, 4*size)

    return CoeffMatrix(Ta, A, B, Cm, dC, D1, Vxx, Vyy, Vzz, dVzz, d2Vzz, Vxy, Vxz, dVxz, Vyz, dVyz)
end

"""
    assemble_direct_matrices(cof::CoeffMatrix, D, D2, be, omega, R)

Assemble the coefficient matrices for the direct quadratic eigenvalue problem:
    L2 * α^2 + L1 * α + L0 = 0

where α is the complex streamwise wavenumber.
"""
function assemble_direct_matrices(cof::CoeffMatrix, D, D2, be, omega, R)
    L0 = cof.D1 + im * R * be * cof.B - im * omega * cof.Ta -
         be^2 * R^2 * cof.Vyy + (cof.C .+ im * be * R * cof.Vyz) * kron(I(4), D) +
         cof.Vzz * kron(I(4), D2)
    L1 = im * cof.A - be * R * cof.Vxy + im * cof.Vxz * kron(I(4), D)
    L2 = -cof.Vxx
    return L0, L1, L2
end

"""
    assemble_adjoint_matrices(cof::CoeffMatrix, D, D2, be, omega, R)

Assemble the coefficient matrices for the adjoint quadratic eigenvalue problem.
"""
function assemble_adjoint_matrices(cof::CoeffMatrix, D, D2, be, omega, R)
    A0 = transpose(cof.D1) + (im * be * R * transpose(cof.B)) -
         (im * omega * transpose(cof.Ta)) - (be^2 * R^2 * transpose(cof.Vyy)) -
         transpose(cof.dC) - (im * be * R * transpose(cof.dVyz)) +
         transpose(cof.d2Vzz) -
         (transpose(cof.C) + im * be * R * transpose(cof.Vyz) - 2 * transpose(cof.dVzz)) * kron(I(4), D) +
         transpose(cof.Vzz) * kron(I(4), D2)
    A1 = (im * transpose(cof.A)) - (be * R * transpose(cof.Vxy)) -
         (im * transpose(cof.dVxz)) - (im * transpose(cof.Vxz)) * kron(I(4), D)
    A2 = -transpose(cof.Vxx)
    return A0, A1, A2
end
