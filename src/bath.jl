export

    bath,
    bath2,
    bath_exact,
    bath2_exact,
    part_func,
    expv
    
    

"""
$(TYPEDSIGNATURES)
"""

function expv(τ::Number, ham::Union{BoseHubbard{T},RBoseHubbard{T}}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(
        ham.H, τ, dense(v, ham.basis), ishermitian=true, tol=1E-8
    )
    @assert info.converged == 1
    U_dket
end
"""
$(TYPEDSIGNATURES)
"""

"""
$(TYPEDSIGNATURES)
"""

function bath(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

#therm_ket = expv(-beta, H[2], state)
    evol_bra = expv( s, H[2], state)
    half_bra = dense(create( State(evol_bra, H[2].basis), j), H[1].basis)
    


    
    evol_ket = expv(τ, H[2], state)
    half_ket = expv(-(τ-s), H[1], create(State(evol_ket, H[2].basis),i))
    
    
    dot(half_bra,half_ket)
   
end

function bath_exact(
    time1::T, time2::T, H::Vector{RBoseHubbard{S}}, i::Int, j::Int, state::State, rho; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

    #therm_ket = expv(-beta, H[2], state)
    evol_bra = expv( s, H[2], state)
    half_bra = dense(create( State(evol_bra, H[2].basis), j), H[1].basis)
    
    
    ket = Matrix(rho)*dense(state, H[2].basis)
    evol_ket = expv(τ, H[2], State(ket, H[2].basis))
    half_ket = expv(-(τ-s), H[1], create(State(evol_ket, H[2].basis),i))

    
    
    
    dot(half_bra,half_ket)
end


function bath2(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

#therm_ket = expv(-beta, H[2], state)
    evol_bra = expv( s, H[2], state)
    print(evol_bra)
    half_bra = dense(destroy( State(evol_bra, H[2].basis), j), H[3].basis)


    
    evol_ket = expv(τ, H[2], state)
    half_ket = expv(-(τ-s), H[3], destroy(State(evol_ket, H[2].basis),i))




dot(half_bra,half_ket)
end   

function bath2_exact(
    time1::T, time2::T, H::Vector{RBoseHubbard{S}}, i::Int, j::Int, state::State, rho; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

#therm_ket = expv(-beta, H[2], state)
    evol_bra = expv( s, H[2], state)
    half_bra = dense(destroy( State(evol_bra, H[2].basis), j), H[2].basis)


    ket = Matrix(rho)*dense(state, H[2].basis)
    evol_ket = expv(τ, H[2], State(ket, H[2].basis))
    half_ket = expv(-(τ-s), H[2], destroy(State(evol_ket, H[2].basis),i))




dot(half_bra,half_ket)

end

function part_func(
    beta::T, H::Vector{BoseHubbard{S}}; kwargs=()
    ) where{S, T <:Real}                                        
    
    eigenvals, _ = eigen(Matrix(H[2].H))
    part_sum = 0
    for (_, energy) in enumerate(eigenvals)
        part_sum = part_sum + exp(-beta*energy) #compute gibbs factor for each eigenstate
    end
    
    part_sum
end 


"""

function thermal_corr(beta::T, H::Vector{BoseHubbard{S}}, time::T; kwargs=()
) where{S, T <:Real}
    

    trsum1 = 0
    trsum2 = 0
    for i in 1:length(eigenvals)
        
        gamma1 =bath(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis), beta )
        gamma2 = bath2(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis), beta )
       
        trsum1 = trsum1 + gamma1
        trsum2 = trsum2 + gamma2
        
    end
    
    corr_arr = [trsum1/ partition_function, trsum2/partition_function]
return corr_arr
end
"""

"""
function bath(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

    
    evol_ket = expv((τ), H[2],state)
    half_ket = expv(-(τ-s), H[1], create(State(evol_ket,H[2].basis), i))
    
    evol_bra = expv( s, H[2],state)

    half_bra = dense(create(State(evol_bra,H[2].basis),j), H[1].basis)    

    
    
    
    dot(half_bra,half_ket)
end


function bath2(
    time1::T, time2::T, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State; kwargs=()
) where {S, T <: Real}
    τ = -1im * time1
    s = -1im * time2

    evol_ket = expv((τ),H[2],state)
    
    a_i_ket= destroy(State(evol_ket,H[2].basis),i)
                
    half_ket = expv(-(τ-s),H[3],a_i_ket)

    #U_ai_ket = dense(State(expv(-(τ-s), H[2], create(State(state_temp,H[3].basis), j)),H[2].basis), H[2].basis)
   
    evol_bra = expv( (s) ,H[2],state)
    half_bra = dense(destroy(State(evol_bra,H[2].basis),j),  H[3].basis)
    
    
    
    dot(half_bra, half_ket)
end
"""



