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
        println(ket,Cₙ[k])
        vecs[k] = ( ket[site] > 0 && Cₙ[k] > 0.001) ? destroy(ket, site) : 0
    end
    print(vecs)
    Ψout = zeros(length(basis_des))
    pos = findall( !iszero, vecs)
    for i in pos
        index = find_index( vecs[i], basis_des)
        Ψout[index] = Cₙ[i]
    end

    return Ψout
end

export super_destroy

