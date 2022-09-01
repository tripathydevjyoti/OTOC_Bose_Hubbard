


using LinearAlgebra
using JeszenszkiBasis
using GraphPlot
using ExponentialAction
using Plots
using KrylovKit



#include("time_evolution.jl")
include("create_destroy_index.jl")
#include("superposition.jl")
include("sparse_hamiltonian.jl")
#include("big_create.jl")
#include("big_destroy.jl")
include("super destory,create.jl")

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
    #res = index
    track = complex(track)
    
    
    time = -1im*time
    
    Ψ₂, info = exponentiate(H2, time, track, ishermitian = true)
    @assert info.converged == 1


    #cₙ,ψₙ = superposition(Ψ₂,basis2)
    #res = ψₙ
  
    #Ψ₃ = big_destroy(cₙ,ψₙ,basis3,site2)

    Ψ₃ = super_destroy(Ψ₂, basis2, basis3, site2)




    Ψ₄, info = exponentiate(H3, -time, Ψ₃, ishermitian = true, kwargs...)
    @assert info.converged == 1

    #cₙ,ψₙ = superposition(Ψ₄,basis3)
    #Ψ₅ = big_create(cₙ,ψₙ,basis2,site1)

    Ψ₅ = super_create(Ψ₄, basis3, basis2, site1)
    


  
    Ψ₆, info = exponentiate(H2, time, Ψ₅, ishermitian = true, kwargs...)
    @assert info.converged == 1

  

    Ψ₇ = super_create(Ψ₆, basis2, basis1, site2)
    

    
    

    Ψ_fin, info = exponentiate(H1, -time, Ψ₇, ishermitian = true, kwargs...)
    @assert info.converged == 1

    track = zeros(length(basis1))
    index = find_index(Ψ₀,basis1)
    track[index] = 1
    

   # return abs(Ψ_fin[index])
    return abs(dot(track, Ψ_fin)) 
    #return Ψ_fin
 
    #return  res
       



export OTOC_lattice
OTOC_lattice(6,1,1.0)
t_vals = range(0,0.5,40)
vals = OTOC_lattice.(5,1,t_vals)
vals =[]
for t in t_vals
    push!(vals, OTOC_lattice(5,1,t))
end    

scatter(t_vals,vals)


print(vals)
end

using BoseHubbardDiagonalize
using JeszenszkiBasis
using LightGraphs #might be ueful later for finding nearest neighbour on hex lattice
using GraphPlot
using SparseArrays: sparse
using DocStringExtensions

function sparse_hamiltonian_test(basis::AbstractSzbasis,L:: Int)
    rows=Int64[]
    cols=Int64[]
    elements=Float64[]
    U=16
    t=4
    

    for (i,bra) in enumerate(basis)
        #diagonal entries
        Usum=0
        for j in 1:basis.K
            Usum += bra[j]*(bra[j]-1)
        end
        push!(rows,i)
        push!(cols,i)
        push!(elements,U*Usum/2)
        println(U*Usum/2)
   
        #off-diagonal entries

        for j in 1:L
            #j_next = mod(j+1, 1:basis.K)
            j_next = j+1
            if j == L 
                j_next = 1
            end
            
              

            for (site1,site2) in [(j,j_next),(j_next,j)]
                if bra[site1]>0
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis
                        push!(rows, i)
                        push!(cols, serial_num(basis,ket))
                        push!(elements, -t*sqrt(bra[site1])*sqrt(bra[site2]+1))
                    end
                end
            end
        end
    end

    sparse(rows,cols,elements,length(basis),length(basis))

end



  







   



 

        
     





basis = Szbasis(2,3)
for i in basis
    println(i)
end
H = Matrix(sparse_hamiltonian_test(basis, 2))
Hexp = exp(+1im*H)
v = *(Hexp,[-0.16074307339206692 - 0.5501089509250516im,
0.5293984495100408 + 0.31847790485748123im,
0.1403431787957471 - 0.0327890897719198im,0])

abs(dot([0,0,1,0],v))
eigvecs(H)
eigvals(H)

L =2
N =3
t = 1.0

ψ₀ = [1,2]


basis1 = Szbasis(L,N)
basis2 = Szbasis(L,N-1)
basis3 = Szbasis(L,N-2)


H1 = sparse_hamiltonian(basis1, L)
H2 = sparse_hamiltonian(basis2, L)
H3 = sparse_hamiltonian(basis3, L)

destroy([1,2],2)
find_index(destroy([1,2],2),basis2)
m2 = exp(-1im*Matrix(H2))
v = m2*[0,1,0]

m3 = exp(+1im*Matrix(H3))
v = m3*[-0.5129516689952816 + 0.07543727000378325im,
-0.4460381033994975 - 0.5132539502570286im]

super_create(v, basis3, basis2, 2)