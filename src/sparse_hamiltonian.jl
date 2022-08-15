using BoseHubbardDiagonalize
using JeszenszkiBasis
using LightGraphs #might be ueful later for finding nearest neighbour on hex lattice
using GraphPlot
using SparseArrays: sparse
using DocStringExtensions

function sparse_hamiltonian(basis::AbstractSzbasis)
    rows=Int64[]
    cols=Int64[]
    elements=Float64[]

    for (i,bra) in enumerate(basis)
        #diagonal entries
        Usum=0
        for j in 1:basis.K
            Usum += bra[j]*(bra[j]-1)
        end
        push!(rows,i)
        push!(cols,i)
        push!(elements,U*Usum/2)

        #off-diagonal entries

        for j in 1:L
            j_next = mod(j+1, 1:basis.K)

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


function time_evolution(
    Eₙ::Array{T, 1},
    ψₙ::Array{T, 2},
    ψ₀::Array{Complex{T}, 1},
    t::T
) where T <: Real
    cₙ = dot.(conj.(ψ₀), ψₙ)
    eₙ = exp.(-1im .* Eₙ .* t)
    sum(eₙ .* cₙ .* ψₙ, dims=2)
end


function superposition(
    Ψ::Array{Complex{T}, 1},
    basis::Array{Int64, 1}
) where T <: Real
    Ψₙ = Int64[]
    β = Complex[]
    
    for i in 1:length(Ψ)
        
        Cₙ = abs(Ψ[i])
        if Cₙ>0.001
            push!(Ψₙ,basis[i])
            push!(β,Cₙ)
        end
        
    end
    return β,Ψₙ
end



function create(
    ket::Array{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] += 1
    return nket
end

function destroy(
    ket::Array{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] -= 1
    return nket
end

function find_index(
    inp_state::Array{Int64,1},
    basis::Array{int64,1}
    )
    index=0
    
    for i in 1:length(basis)
        if inp_state == basis[i]
            index =i
            break
        end

    end
    

    return index
end

#find_index([5,1,0,0,0,0],basis1)
   

function big_destroy(
    Cₙ::Array{Complex{T},1},
    Ψin::Array{Int64,2},
    basis::Array{Int64,1},
    site::Int64
) where T <:Real
    
    Ψout = zeros(length(basis))

    for (i,ϕ) in enumerate(Ψin)
        ϕ = destroy( ϕ , site)
        index = find_index(ϕ,basis)
        Ψout[index] = Cₙ[i]
    end
    return Ψout
end    

function big_create(
    Cₙ::Array{Complex{T},1},
    Ψin::Array{Int64,2},
    basis::Array{Int64,1},
    site::Int64
) where T <:Real

    Ψout = zeros(length(basis))

    for (i,ϕ) in enumerate(Ψin)
        ϕ = create( ϕ , site)
        index = find_index(ϕ,basis)
        Ψout[index] = Cₙ[i]
    end
    return Ψout
end        

        
     


function OTOC(
    site1::Int64,
    site2::Int64,
    time::Real)
    
     Ψ₀=[1,1,1,1,1,1]

     Ψ₁ = destroy(Ψ₀,site1)

    index::Int64 = find_index(Ψ₁,basis2)
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
