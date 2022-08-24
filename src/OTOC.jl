module OTOC_example


using LinearAlgebra
using JeszenszkiBasis
using GraphPlot
using ExponentialAction
using Plots
using KrylovKit



include("time_evolution.jl")
include("create_destroy_index.jl")
include("superposition.jl")
include("sparse_hamiltonian.jl")
include("big_create.jl")
include("big_destroy.jl")

export sparse_hamiltonian
export time_evolution
export find_index

basis1 = Szbasis(6,6)

basis2 = Szbasis(6,5)

basis3 = Szbasis(6,4)


H1 = sparse_hamiltonian(basis1)
      
H2 = sparse_hamiltonian(basis2)
H3 = sparse_hamiltonian(basis3)





function OTOC_lattice(
    site1::Int64,
    site2::Int64,
    time::Real)

    #E1 = E1.values

    #E2 = E2.values

    #E3 = E3.values




    kwargs = ()
    Ψ₀=[1,1,1,1,1,1]

    Ψ₁ = destroy(Ψ₀,site1)

    index = find_index(Ψ₁,basis2)
    track = zeros(length(basis2))
    track[index]=1
    res = index
    #track = complex(track)
    
    #Ψ₂ = time_evolution(E2.values,V2,track,time)
    time = -1im*time
    #Ψ₂ = expv(time, Matrix(H2) , track)
    Ψ₂, info = exponentiate(H2, time, track, ishermitian = true)
    @assert info.converged == 1


    cₙ,ψₙ = superposition(Ψ₂,basis2)
    res = ψₙ
  
    Ψ₃ = big_destroy(cₙ,ψₙ,basis3,site2)


    #Ψ₃ = complex(Ψ₃)
    #Ψ₄ = time_evolution(E3.values,V3,Ψ₃,-time)
    #Ψ₄ = expv(-time, Matrix(H3), Ψ₃)

    Ψ₄, info = exponentiate(H3, -time, Ψ₃, ishermitian = true, kwargs...)
    @assert info.converged == 1

    cₙ,ψₙ = superposition(Ψ₄,basis3)
    Ψ₅ = big_create(cₙ,ψₙ,basis2,site1)
    #Ψ₅ = complex(Ψ₅)


    #Ψ₆ = time_evolution(E2.values,V2,Ψ₅,time)
    #Ψ₆ = expv(time, Matrix(H2), Ψ₅)
    Ψ₆, info = exponentiate(H2, time, Ψ₅, ishermitian = true, kwargs...)
    @assert info.converged == 1

    cₙ,ψₙ = superposition(Ψ₆,basis2) 
    Ψ₇ = big_create(cₙ,ψₙ,basis1,site2)
    #Ψ₇ = complex(Ψ₇)

    #Ψ_fin = time_evolution(E1.values,V1,Ψ₇,-time)
    Ψ_fin = expv(-time, Matrix(H1), Ψ₇)

    
    index = find_index(Ψ₀,basis1)
    

    return abs(Ψ_fin[index])
   
    #return Ψ_fin
 
    #return  res
       

end

export OTOC_lattice
OTOC_lattice(5,1,1.0)
t_vals = range(0,0.8,100)
vals = OTOC_lattice.(6,1,t_vals)
vals =[]
for t in t_vals
    push!(vals, OTOC_lattice(5,1,t))
end    

scatter(t_vals,vals)
end