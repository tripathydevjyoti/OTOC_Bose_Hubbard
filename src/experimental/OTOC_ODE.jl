export
    OTOC_ODE

"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Number, ham::BoseHubbard{T}, state::State) where T
    f(du, u, p, t) = mul!(du, c .* ham.H, u)
    ket = convert.(Complex{T}, dense(state, ham.basis))
    sol = solve(ODEProblem(f, ket, τ), Tsit5())
    sol(τ)
end

"""
$(TYPEDSIGNATURES)
"""
function OTOC_ODE(
    H::Vector{BoseHubbard{S}}, i::Int, j::Int, ket::State, time::Vector{T}
) where {S, T <: Number}
    otoc = Complex{T}[]

    ai_ket = destroy(ket, i)

    for τ ∈ time
        # 1. |x> := V * a_j * U * a_i * ket
        U_ai_ket = expv(τ, -1im, H[2], ai_ket)
        aj_U_ai_ket = destroy(State(U_ai_ket, H[2].basis), j)
        V_aj_U_ai_ket = expv(τ, 1im, H[3], aj_U_ai_ket)

        # 2. |y> := a_i * V * a_j * U * ket
        U_ket = expv(τ, -1im, H[1], ket)
        aj_U_ket = destroy(State(U_ket, H[1].basis), j)
        V_aj_U_ket = expv(τ, 1im, H[2], aj_U_ket)
        ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

        # 3. OTOC: <y|x>
        push!(otoc, dot(ai_V_aj_U_ket, V_aj_U_ai_ket))
    end
    otoc
end
