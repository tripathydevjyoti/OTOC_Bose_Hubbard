export
    OTOC

"""
$(TYPEDSIGNATURES)
"""
function OTOC(
    ham::BoseHubbard{S}, i::Int, j::Int, state::State, time::T; kwargs=()
) where {S, T <: Real}
    τ = -1im * time

    # 1. compute |x> := V * a_j * U * a_i * ket
    ai_ket = dense(destroy(state, i), ham.basis)
    U_ai_ket, info = exponentiate(ham.H, τ, ai_ket, ishermitian=true)
    @assert info.converged == 1

    aj_U_ai_ket = dense(destroy(State(U_ai_ket, ham.basis), j), ham.basis)
    V_aj_U_ai_ket, info = exponentiate(ham.H, -τ, aj_U_ai_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    # 2. compute |y> := a_i * V * a_j * U * ket
    U_ket, info = exponentiate(ham.H, τ, dense(state, ham.basis), ishermitian=true, kwargs...)
    @assert info.converged == 1

    aj_U_ket = dense(destroy(State(U_ket, ham.basis), j), ham.basis)
    V_aj_U_ket, info = exponentiate(ham.H, -τ, aj_U_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, ham.basis), i), ham.basis)

    # 3. compute OTOC: <y|x>
    dot(ai_V_aj_U_ket, V_aj_U_ai_ket)
end

function OTOC(
    ham::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State, time::T; kwargs=()
) where {S, T <: Real}
    τ = -1im * time

    # 1. compute |x> := V * a_j * U * a_i * ket
    ai_ket = dense(destroy(state, i), ham[2].basis)
    U_ai_ket, info = exponentiate(ham[2].H, τ, ai_ket, ishermitian=true)
    @assert info.converged == 1

    aj_U_ai_ket = dense(destroy(State(U_ai_ket, ham[2].basis), j), ham[3].basis)
    V_aj_U_ai_ket, info = exponentiate(ham[3].H, -τ, aj_U_ai_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    # 2. compute |y> := a_i * V * a_j * U * ket
    U_ket, info = exponentiate(ham[1].H, τ, dense(state, ham[1].basis), ishermitian=true, kwargs...)
    @assert info.converged == 1

    aj_U_ket = dense(destroy(State(U_ket, ham[1].basis), j), ham[2].basis)
    V_aj_U_ket, info = exponentiate(ham[2].H, -τ, aj_U_ket, ishermitian=true, kwargs...)
    @assert info.converged == 1

    ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, ham[2].basis), i), ham[3].basis)

    # 3. compute OTOC: <y|x>
    dot(ai_V_aj_U_ket, V_aj_U_ai_ket)
end
