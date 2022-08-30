module OTOC_example



using LinearAlgebra
using JeszenszkiBasis
using GraphPlot
using ExponentialAction
using Plots
using KrylovKit




include("create_destroy_index.jl")
include("sparse_hamiltonian.jl")
include("super destory,create.jl")









function OTOC_lattice(
    site1::Int64,
    site2::Int64,
    time::Real
)

    kwargs = ()
    Ψ₀=[1,1,1,1,1,1]

    Ψ₁ = destroy(Ψ₀,site1)

    index = find_index(Ψ₁,basis2)
    track = zeros(length(basis2))
    track[index]=1
    track = complex(track)
    
    
    time = -1im*time
    
    Ψ₂, info = exponentiate(H2, time, track, ishermitian = true)
    @assert info.converged == 1

    Ψ₃ = super_destroy(Ψ₂, basis2, basis3, site2)

    Ψ₄, info = exponentiate(H3, -time, Ψ₃, ishermitian = true, kwargs...)
    @assert info.converged == 1


    Ψ₅ = super_create(Ψ₄, basis3, basis2, site1)
    
    Ψ₆, info = exponentiate(H2, time, Ψ₅, ishermitian = true, kwargs...)
    @assert info.converged == 1

  

    Ψ₇ = super_create(Ψ₆, basis2, basis1, site2)
    
    Ψ_fin, info = exponentiate(H1, -time, Ψ₇, ishermitian = true, kwargs...)
    @assert info.converged == 1

    track = zeros(length(basis1))
    index = find_index(Ψ₀,basis1)
    track[index] = 1
    

    return abs(Ψ_fin[index])
    #return real(dot(track, Ψ_fin)) 
    
       

end

export OTOC_lattice
end


