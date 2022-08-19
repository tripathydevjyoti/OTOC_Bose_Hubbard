using JeszenszkiBasis

function big_destroy(
    Cₙ::Vector{Complex},
    Ψin ::Vector{AbstractArray{Int64}},
    basis::AbstractSzbasis,
    site::Int64
) 
    
    Ψout = zeros(length(basis))

    for (i,ϕ) in enumerate(Ψin)
        ϕ = destroy( ϕ , site)
        index = find_index(ϕ,basis)
        Ψout[index] = Cₙ[i]
    end
    return Ψout
end