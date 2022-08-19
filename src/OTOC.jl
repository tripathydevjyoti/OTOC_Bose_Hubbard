module OTOC

using LinearAlgebra
using JeszenszkiBasis
using GraphPlot



include("time_evolution.jl")
include("create_destroy_index.jl")
include("superposition.jl")
include("sparse_hamiltonian.jl")
include("big_create.jl")
include("big_destroy.jl")



function OTOC_lattice(
    site1::Int64,
    site2::Int64,
    time::Real)
    
    Ψ₀=[1,1,1,1,1,1]

    Ψ₁ = destroy(Ψ₀,site1)

    index = find_index(Ψ₁,basis2)
    track = zeros(length(basis2))
    track[index]=1

    Ψ₂ = time_evolution(E2,V2,track,time)

    cₙ,ψₙ = superposition(Ψ₂,basis2)
    Ψ₃ = big_destroy(cₙ,ψₙ,basis3,site2)

    Ψ₄ = time_evolution(E3,V3,Ψ₃,-time)

    cₙ,ψₙ = superposition(Ψ₄,basis3)
    Ψ₅ = big_create(cₙ,ψₙ,basis2,site1)

    Ψ₆ = time_evolution(E2,V2,Ψ₅,time)

    cₙ,ψₙ = superposition(Ψ₆,basis2) 
    Ψ₇ = big_create(cₙ,ψₙ,basis1,site2)

    Ψ_fin = time_evolution(E1,V1,Ψ₇,-time)

    
    index = find_index(init_state1,basis1)
    
    
    return abs(Ψ_fin[index])

end

end