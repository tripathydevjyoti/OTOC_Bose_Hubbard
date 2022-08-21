export
    operate,
    create, destroy,
    creation, annihilation

"""
$(TYPEDSIGNATURES)
"""
function operate(ket::Vector{Int}, i::Int, op::Symbol)
    @assert op ∈ (:create, :destroy)
    nket = copy(ket)
    op == :destroy ? nket[i] -= 1 : nket[i] += 1
    nket
end
create(ket::Vector{Int}, i::Int) = operate(ket, i, :create)
destroy(ket::Vector{Int}, i::Int) = operate(ket, i, :destroy)

"""
$(TYPEDSIGNATURES)
"""
function annihilation(::Type{T}, B::Basis, i::Int) where T <: Real
    n = length(B.eig_vecs)
    I, J, V = Int[], Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        if ket[i] > 0
            push!(J, v)
            push!(I, get_index(B, operate(ket, i, :destroy)))
            push!(V, T(ket[i]) |> sqrt)
        end
    end
    sparse(I, J, V, n, n)
end
annihilation(B::Basis, i::Int) = annihilation(Float64, B::Basis, i::Int)

"""
$(TYPEDSIGNATURES)
"""
creation(::Type{T}, B::Basis, i::Int) where {T} = transpose(annihilation(T, B, i))
