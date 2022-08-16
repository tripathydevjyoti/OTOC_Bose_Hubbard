function superposition(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis
)   where T<: Real
    Ψₙ = ComplexF64[]
    #Ψₙ = Vector{ComplexF64,1}
    β = Complex[]
    
    Cₙ =broadcast(abs, Ψ)
    if 
    for i in 1:length(Ψ)
        
        Cₙ = abs(Ψ[i])
        if Cₙ>0.001
            push!(Ψₙ,basis[i])
            push!(β,Ψ[i])
        end
        
    end

    
    return β,Ψₙ 
    
end