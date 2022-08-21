export
    OTOC

"""
$(TYPEDSIGNATURES)
"""
function OTOC(
    ham::S, i::Int, j::Int, state::State, time::T; kwargs=()
) where {S <: BoseHubbard, T <: Real}
    τ = 1im * time

    # 1. compute |x> := V * ap_j * U * a_i * ket
    a_ket = dense(
        State(state.coeff, destroy.(state.eig_vecs, i)), ham.basis
    )
    U_a_ket, info = exponentiate(ham.H, τ, a_ket, ishermitian=true)
    @assert info.converged == 1

    idx = findall(!iszero, U_a_ket)
    ap_U_a_ket = dense(
        State(U_a_ket[idx], create.(ham.basis.eig_vecs[idx], j)), ham.basis
    )
    V_ap_U_a_ket, info = exponentiate(ham.H, -τ, ap_U_a_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    # 2. compute |y> := a_i * V * ap_j * U * ket
    U_ket, info = exponentiate(ham.H, τ, dense(state, ham.basis), ishermitian=true, kwargs...)
    @assert info.converged == 1

    idx = findall(!iszero, U_ket)
    ap_U_ket = dense(
        State(U_ket[idx], create.(ham.basis.eig_vecs[idx], j)), ham.basis
    )
    V_ap_U_ket, info = exponentiate(ham.H, -τ, ap_U_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    idx = findall(!iszero, V_ap_U_ket)
    a_V_ap_U_ket = dense(
        State(V_ap_U_ket[idx], destroy.(ham.basis.eig_vecs[idx], j)), ham.basis
    )

    # 3. <y|x>
    dot(a_V_ap_U_ket, V_ap_U_a_ket)
end
