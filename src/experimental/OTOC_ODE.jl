export
    OTOC_ODE

"""
$(TYPEDSIGNATURES)
"""
function ode_expv(τ::Real, c::Number, ham::BoseHubbard{T}, state::State) where T
    ket = dense(state, ham.basis)
    z = zero(eltype(τ))

    if τ ≈ z return ket end
    A(du, u, p, t) = mul!(du, c .* ham.H, u)

    ODE = ODEProblem(A, ket, (z, τ))
    alg = RK4()
    #alg = ETDRK4(krylov=true, m=30)

    sol = solve(
        ODE, alg, save_everystep=false, reltol=1e-8, abstol=1e-8
    )

    Uket = sol[end]
    Uket ./ norm(Uket)
end

"""
$(TYPEDSIGNATURES)
"""
function OTOC_ODE(
    times::Vector{T}, H::Vector{BoseHubbard{S}}, i::Int, j::Int, ket::State
) where {S, T <: Number}
    otoc = Complex{T}[]

    ai_ket = destroy(ket, i)
    for τ ∈ times
        # 1. compute |x> := V * a_j * U * a_i * ket
        U_ai_ket = ode_expv(τ, -1im, H[2], ai_ket)
        V_aj_U_ai_ket = ode_expv(τ, 1im, H[3], destroy(State(U_ai_ket, H[2].basis), j))

        # 2. compute |y> := a_i * V * a_j * U * ket
        U_ket = ode_expv(τ, -1im, H[1], ket)
        V_aj_U_ket = ode_expv(τ, 1im, H[2], destroy(State(U_ket, H[1].basis), j))
        ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

        # 3. compute OTOC: <y|x>
        push!(otoc, dot(ai_V_aj_U_ket, V_aj_U_ai_ket))
    end
    otoc
end
