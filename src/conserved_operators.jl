export
    occupation

function occupation(::Type{T}, B::Basis) where T <: Real
    n = length(B.eig_vecs)
    I, V = Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        push!(I, v)
        push!(V, sum(ket))
    end
    sparse(I, I, V, n, n)
end

function occupation(::Type{T}, B::Basis, i::Int) where T
    n = length(B.eig_vecs)
    I, V = Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        push!(I, v)
        push!(V, ket[i])
    end
    sparse(I, I, V, n, n)
end

occupation(B::Basis) = occupation(Float64, B::Basis)
occupation(B::Basis, i::Int) = occupation(Float64, B::Basis, i)
