using JeszenszkiBasis

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