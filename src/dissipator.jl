export
    diss_one_cre,
    diss_one_des,
    diss_two_cre,
    diss_two_des

function expv(τ::Number, ham::BoseHubbard{T}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(ham.H, τ, dense(v, ham.basis), ishermitian = true, tol = 1E-8)
    @assert info.converged == 1
    U_dket
end

function diss_one_cre(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs...
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime

    U_ket = expv(τ, H[2], ket)
    ai_dag_U_ket = create(State(U_ket,H[2].basis),i)
    Udag_ai_dag_U_ket = expv(-τ, H[2], ai_dag_U_ket )
    mid_ket = State(rho_t * Udag_ai_dag_U_ket, H[2].basis)

    U_mid_ket = expv(τ-s, H[2], mid_ket)
    Udag_aj_U_mid_ket = expv(s-τ, H[2], destroy(State(U_mid_ket, H[2].basis),j))

    dot(dense(bra, H[2].basis), Udag_aj_U_mid_ket)
end

function diss_one_des(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs...
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime


    U_ket = expv(τ, H[2], ket)
    ai_U_ket = destroy(State(U_ket,H[2].basis),i)
    Udag_ai_U_ket = expv(-τ, H[2], ai_U_ket )
    mid_ket = State(rho_t * Udag_ai_U_ket, H[2].basis)

    U_mid_ket = expv(τ-s, H[2], mid_ket)
    Udag_aj_dag_U_mid_ket = expv(s-τ, H[2], create(State(U_mid_ket, H[2].basis),j) )

    dot(dense(bra, H[2].basis), Udag_aj_dag_U_mid_ket)
end

function diss_two_cre(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs...
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime

    rho_ket = State(rho_t * dense(ket,H[2].basis), H[2].basis)
    U_rho_ket = expv(τ-s, H[2], rho_ket)
    U_dag_aj_U_rho_ket = expv(-s, H[2], destroy(State(U_rho_ket, H[2].basis),j))
    mix_U_state = U_dag_aj_U_rho_ket
    Udag_ai_dag_mix_U_state = expv(-τ, H[2], create(State(mix_U_state, H[2].basis), i))

    dot(dense(bra,H[2].basis), Udag_ai_dag_mix_U_state)
end

function diss_two_des(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs...
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime

    rho_ket = State(rho_t * dense(ket,H[2].basis), H[2].basis)
    U_rho_ket = expv(τ-s, H[2], rho_ket)
    U_dag_aj_dag_U_rho_ket = expv(-s, H[2], create(State(U_rho_ket, H[2].basis),j))
    mix_U_state = U_dag_aj_dag_U_rho_ket
    Udag_ai_mix_U_state = expv(-τ, H[2], destroy(State(mix_U_state, H[2].basis),i))

    dot(dense(bra,H[2].basis), Udag_ai_mix_U_state)
end
