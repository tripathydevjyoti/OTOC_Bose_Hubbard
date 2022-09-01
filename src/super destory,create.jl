using JeszenszkiBasis
function super_destroy(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis,
    basis_des ::AbstractSzbasis,
    site ::Int
) where T<:Real 

    Cₙ = abs.(Ψ)
    n = length(Cₙ)
    vecs = Vector(undef, n)
    for k in 1:n
        ket = basis[k]
        vecs[k] = ( ket[site] > 0 && Cₙ[k] > 0.001) ? destroy(ket, site) : 0
    end
    
    Ψout = complex(zeros(length(basis_des)))
    pos = findall( !iszero, vecs)
    for i in pos
        index = find_index( vecs[i], basis_des)
        Ψout[index] = Ψ[i]
    end

    return Ψout
end
export super_destroy




function super_create(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis,
    basis_cre ::AbstractSzbasis,
    site ::Int
) where T<:Real 

    Cₙ = abs.(Ψ)
    n = length(Cₙ)
    vecs = Vector(undef, n)
    for k in 1:n
        ket = basis[k]
        
        vecs[k] = Cₙ[k] > 0.001 ? create(ket, site) : 0
    end
    
    Ψout = complex(zeros(length(basis_cre)))
    pos = findall( !iszero, vecs)
    
    for i in pos
        index = find_index( vecs[i], basis_cre)
        Ψout[index] = Ψ[i]
    end

    return Ψout
end
export super_create


