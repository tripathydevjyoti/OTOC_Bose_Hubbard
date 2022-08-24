using JeszenszkiBasis



function big_destroy(
    Cₙ::Vector{Complex},
   
    Ψin ::Vector{Vector{Int64}},
    basis::AbstractSzbasis,
    site::Int64
) 
    
    Ψout = zeros(length(basis))
    Ψout = complex(Ψout)
    

    for (i,ϕ) in enumerate(Ψin)
        
        pre_sum = sum(ϕ)
        ϕ = destroy( ϕ , site)
        
       
        if sum(ϕ) == pre_sum - 1
            index = find_index(ϕ,basis)
            Ψout[index] = Cₙ[i]
            print(Cₙ[i])

        end

        
    end
    return Ψout
end

"""
basis = Szbasis(3,2)
Ψ = [1/5, 1/5, 1/5, 1/5, 1/5]
Ψ = complex(Ψ)

coeff_list, states_list = superposition(Ψ, basis)
coeff_list
states_list
for (i,x) in enumerate(states_list)
    println(i,x)
end    
basis1 =Szbasis(3,1)
    
big_destroy(coeff_list,states_list,basis1,1)

"""
    
