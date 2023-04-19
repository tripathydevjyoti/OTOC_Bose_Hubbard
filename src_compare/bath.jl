export

    bath,
    bath2

"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Number, ham::BoseHubbard{T}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(
        ham.H, τ, dense(v, ham.basis), ishermitian=true, tol=1E-8
    )
    @assert info.converged == 1
    U_dket
end

"""
$(TYPEDSIGNATURES)
"""


function bath(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

    
    evol_ket = expv((τ-s), H[2],state)
    half_ket = expv(s, H[3], destroy(State(evol_ket,H[2].basis), i))
    
    evol_bra = expv(τ, H[2],state)

    half_bra = dense(destroy(State(evol_bra,H[2].basis),j), H[3].basis)    

    
    
    
    dot(half_bra,half_ket)
end


function bath2(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

    evol_ket = expv((τ-s),H[2],state)
    
    adag_i_ket= create(State(evol_ket,H[2].basis),i)
                
    state_temp = expv(s,H[1],adag_i_ket)

    U_ai_ket = dense(State(expv(-τ, H[2], destroy(State(state_temp,H[1].basis), j)),H[2].basis), H[2].basis)
   
   
    
    
    
    dot(dense(state,H[2].basis), U_ai_ket)
end




