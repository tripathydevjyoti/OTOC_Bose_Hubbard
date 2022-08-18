using Test
using LinearAlgebra

function time_evolution(
    Eₙ::Array{T, 1},
    #ψₙ::Array{T, 2},
    ψₙ::Vector{AbstractArray{Int64}},
    ψ₀::Array{Complex{T}, 1},
    t::T
) where T <: Real
    cₙ = dot.(conj.(ψ₀), ψₙ)
    eₙ = exp.(-1im .* Eₙ .* t)
    sum(eₙ .* cₙ .* ψₙ, dims=2)
end

e = [0.5,-0.5]


v = [[1,0],[0,1]]
    
Ψ₀ = [1/sqrt(2), 1/sqrt(2)]
Ψ₀ = complex(Ψ₀)




time_evolution(e, v, Ψ₀, 1.0)

@testset "Checking if time_eolution is working fine" begin
    e = [0.5,-0.5]


    v = [[1,0],[0,1]]
    
    Ψ₀ = [1/sqrt(2), 1/sqrt(2)]
    Ψ₀ = complex(Ψ₀)

    time_evolution(e, v, Ψ₀, 1.0)
    
end