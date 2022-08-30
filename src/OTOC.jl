export
    OTOC

"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Number, ham::BoseHubbard{T}, v::State) where T
    U_dket, info = exponentiate(ham.H, τ, dense(v, ham.basis), ishermitian=true)
    @assert info.converged == 1
    U_dket
end

"""
$(TYPEDSIGNATURES)
"""
function OTOC(
    H::BoseHubbard{S}, i::Int, j::Int, state::State, time::T; kwargs=()
) where {S, T <: Real}
    τ = -1im * time

    # 1. compute |x> := V * a_j * U * a_i * ket
    U_ai_ket = expv(τ, H, destroy(state, i))
    V_aj_U_ai_ket = expv(-τ, H, destroy(State(U_ai_ket, H.basis), j))

    # 2. compute |y> := a_i * V * a_j * U * ket
    U_ket = expv(τ, H, state)
    V_aj_U_ket = expv(-τ, H, destroy(State(U_ket, H.basis), j))
    ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H.basis), i), H.basis)

    # 3. compute OTOC: <y|x>
    dot(ai_V_aj_U_ket, V_aj_U_ai_ket)
end

function OTOC(
    H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State, time::T; kwargs=()
) where {S, T <: Real}
    τ = -1im * time

    # 1. compute |x> := V * a_j * U * a_i * ket
    U_ai_ket = expv(τ, H[2], destroy(state, i))
    V_aj_U_ai_ket = expv(-τ, H[3], destroy(State(U_ai_ket, H[2].basis), j))

    # 2. compute |y> := a_i * V * a_j * U * ket
    U_ket = expv(τ, H[1], state)
    V_aj_U_ket = expv(-τ, H[2], destroy(State(U_ket, H[1].basis), j))
    ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

    # 3. compute OTOC: <y|x>
    dot(ai_V_aj_U_ket, V_aj_U_ai_ket)
end
