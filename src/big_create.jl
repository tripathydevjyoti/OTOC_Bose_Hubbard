using JeszenszkiBasis

function big_create(
    Cₙ::Vector{Complex},
    
    Ψin::Vector{Vector{Int64}},
    basis::AbstractSzbasis,
    site::Int64
) 

    Ψout = zeros(length(basis))
    Ψout = complex(Ψout)


    for (i,ϕ) in enumerate(Ψin)
        
        ϕ = create( ϕ , site)
        
        index = find_index(ϕ,basis)
        
        Ψout[index] = Cₙ[i]
    end
    return Ψout
end    
"""
basis = Szbasis(3,2)
Ψ = [1/5, 1/5, 1/5]
Ψ = complex(Ψ)
basis1 =Szbasis(3,1)
coeff_list, states_list = superposition(Ψ, basis1)
coeff_list
states_list
for (i,x) in enumerate(states_list)
    println(i,x)
end    

    
big_create(coeff_list,states_list,basis,1)
"""