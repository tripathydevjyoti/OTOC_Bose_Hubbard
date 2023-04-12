export
    diss_one,
    diss_two



function expv(τ::Number, ham::BoseHubbard{T}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(
    ham.H, τ, dense(v, ham.basis), ishermitian=true, tol=1E-8
    )
    @assert info.converged == 1
    U_dket
end    
   

function diss_one(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs=()
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime

    ket_state = State(U_ket,H[2].basis)
    U_ket = expv(τ, H[2], ket_state)
    ai_dag_U_ket = create(State(U_ket,H[2].basis),i)
    Udag_ai_dag_U_ket = expv(-τ, H[1], ai_dag_U_ket )
    mid_ket = rho_t * Udag_ai_dag_U_ket

    U_mid_ket = expv(τ-s, H[1], mid_ket)
    Udag_aj_U_mid_ket = expv(s-τ, H[2], destroy(State(U_mid_ket, H[1].basis),j) )

    dot(bra, Udag_aj_U_mid_ket)

end 

function diss_two(
    time::T, stime::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, bra::State, ket::State, rho_t; kwargs=()
) where {S, T <: Real}
    τ = -1im * time
    s = -1im * stime


    rho_ket = rho_t * ket
    U_rho_ket = expv(τ-s, H[2], rho_ket)
    U_dag_aj_U_rho_ket =expv(s, H[3], destroy(State(U_rho_ket, H[2].basis),j) )
    mix_U_state = U_dag_aj_U_rho_ket
    Udag_ai_dag_mix_U_state = expv(-τ, H[2], create(State(mix_U_state, H[3].basis),i) )


    dot(bra, Udag_ai_dag_mix_U_state)

end    






 