using Test
using LinearAlgebra

function time_evolution(
    Eₙ::Vector{Float64},
    #ψₙ::Array{T, 2},
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

 

@testset "Checking if time_eolution is working fine" begin
    e = [0.5,-0.5]

    c =1/sqrt(2)
    v = [[1,0],[0,1]]
    
    Ψ₀ = [1/sqrt(2), 1/sqrt(2)]
    Ψ₀ = complex(Ψ₀)

    @test time_evolution(e, v, Ψ₀, 1.0) == [ c*exp(-0.5im), c*exp(+0.5im)]
    
end