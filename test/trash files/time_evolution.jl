using Test
using LinearAlgebra

function time_evolution(
    Eₙ::Vector{Float64},
    ψₙ::Vector{Vector{Int64}},
    ψ₀::Array{Complex{T}, 1},
    t::T
) where T <: Real
     
    function inner_product(
        ψ ::Vector{Int64}
    ) 
        δ = dot(Ψ₀,ψ)
    end    

    
    cₙ = inner_product.(v)
    eₙ = exp.(-1im .* e .*t)
    sum(eₙ .* cₙ .* v)
end