export
    OTOC_bose_bubbard

"""
$(TYPEDSIGNATURES)
"""
function OTOC_bose_bubbard(
    H::S, i::Int, j::Int, state::State{S}, time::T; kwargs
) where {S, T <: Real}
    τ = 1im * time
    args = (kwargs..., ishermitian=true)

    # 1. compute |x> := V * ap_j * U * a_i * ket
    a_ket = State(state.coeff, destroy.(state.eig_vecs, i))
    U_a_ket, infoU = exponentiate(H, τ, dense(a_ket, B.basis), args...)
    @assert infoU.converged == 1

    idx = findall(!iszero, U_a_ket)
    ap_U_a_ket = State(U_a_ket[idx], create.(B.basis[idx], j))
    V_ap_U_a_ket, infoV = exponentiate(H, -τ, dense(ap_U_a_ket, B.basis), args...)
    @assert infoV.converged == 1

    # 2. compute |y> := a_i * V * ap_j * U * ket
    U_ket, infoU = exponentiate(H, τ, dense(state), args...)
    @assert infoV.converged == 1

    idx = findall(!iszero, U_ket)
    ap_U_ket = State(U_ket[idx], create.(B.basis[idx], j))
    V_ap_U_ket, infoV = exponentiate(H, -τ, dense(ap_U_ket, B.basis), args...)
    @assert infoV.converged == 1

    idx = findall(!iszero, V_ap_U_ket)
    a_V_ap_U_ket = dense(State(V_ap_U_ket[idx], destroy.(B.basis[idx], j)))

    # 3. <y|x>
    dot(a_V_ap_U_ket, V_ap_U_a_ket)
end
