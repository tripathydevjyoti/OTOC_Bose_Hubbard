
"""
NOT COMPLETE YET


"""

using Test
using JeszenszkiBasis


function find_index(
    inp_state::Array{Int64,1},
    basis::AbstractSzbasis
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


function big_destroy(
    Cₙ::Array{Complex{T},1},
    Ψin::Array{Int64,1},
    basis::AbstractSzbasis,
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




@testset "Checking if Big Create and Destroy are working fine" begin

    c = 1/sqrt(2)
    coeff_list = Float64[]
    push!(coeff_list, c)
    push!(coeff_list, c)

    states_list = Int64[]
    push!(states_list, [2,0,0])
    push!(states_list, [1,1,0])

    basis = Szbasis(3,2)

    print( big_destroy(coeff_list,states_list,basis,1))
    @test big_destroy(coeff_list,states_list,basis,1) == [c,c,0]
end

c = 1/sqrt(2)
coeff_list = Float64[]
push!(coeff_list, c)
push!(coeff_list, c)

basis = Szbasis(3,2)
states_list = Array{Int64,2}
typeof(states_list)
states_list=[ [2,0,0], [1,1,0] ]

print( big_destroy(coeff_list,states_list,basis,1))
