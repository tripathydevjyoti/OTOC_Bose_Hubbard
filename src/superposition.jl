
using Test
using JeszenszkiBasis

function superposition(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis
)   where T<: Real
    

    Ψₙ = Vector{Int64}[]
    β = Complex[]

    Cₙ =broadcast(abs, Ψ)
    
    for(i,x) in enumerate(Cₙ)
        
        
        if x > 0.001
            push!(Ψₙ,basis[i])
            push!(β,Ψ[i])
        end
        
    end
    
    return β,Ψₙ 
    
end

vec = [1,0,0,1]
K= findall(x->x>0.001, vec)

function super_destroy(
    Ψ ::Vector{Complex{T}}
    basis ::AbstractSzbasis
    site ::Int
)   where T<: Real

    Cₙ = abs.(Ψ)
    vecs = Vector(undef, length(Cₙ))


    
   