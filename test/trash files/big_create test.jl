using Test
using JeszenszkiBasis

function create(
    ket::AbstractArray{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] += 1
    return nket
end

function find_index(
    inp_state::Array{Int64,1},
    basis::Array{Int64,1}
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

function superposition(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis
)   where T<: Real
    
    Ψₙ = AbstractArray{Int64}[]
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


function big_create(
    Cₙ::Vector{Complex},
    Ψin::Vector{AbstractArray{Int64}},
    basis::AbstractSzbasis,
    site::Int64
) 

    Ψout = zeros(length(basis))

    for (i,ϕ) in enumerate(Ψin)
        ϕ = create( ϕ , site)
        index = find_index(ϕ,basis)
        Ψout[index] = Cₙ[i]
    end
    return Ψout
end

@testset "Checking if big_create is working fine" begin
    basis = Szbasis(3,1)
    Ψ = [1/5, 4/5, 0]
    Ψ = complex(Ψ)

    coeff_list, states_list = superposition(Ψ, basis)

    basis1 =Szbasis(3,2)
    
    @test big_create(coeff_list,states_list,basis1,1)   == [0.2, 0.8, 0, 0, 0, 0]


    
end
