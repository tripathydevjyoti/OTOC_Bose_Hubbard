export
    operate,
    destroy_and_create,
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

"""
$(TYPEDSIGNATURES)
"""
function destroy_and_create(ket::Vector{Int}, i::Int, j::Int)
    operate(operate(ket, i, :destroy), j, :create)
end

#=
"""
$(TYPEDSIGNATURES)
"""
function annihilation(::Type{T}, B::Basis, i::Int, ::Val{:sparse}) where T <: Real
    n = length(B.eig_vecs)
    I, J, V = Int[], Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        if ket[i] > 0
            push!(J, v)
            tg = tag(operate(ket, i, :destroy))
            push!(I, searchsortedfirst(B.tags, tg))
            push!(V, T(ket[i]) |> sqrt)
        end
    end
    sparse(I, J, V, n, n)
end

function annihilation(::Type{T}, B::Basis, i::Int, ::Val{:dense}) where T <: Real
    n = length(B.eig_vecs)
    a = zeros(T, n, n)
    eig_vecs = enumerate(B.eig_vecs)
    for (v, ket) ∈ eig_vecs
        if ket[i] > 0
            aket = operate(ket, i, :destroy)
            for (w, bra) ∈ eig_vecs
                if bra == aket a[w, v] = T(ket[i]) |> sqrt end
            end
        end
    end
    a
end

annihilation(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = annihilation(T, B, i, Val(s))
annihilation(B::Basis, i::Int, s::Symbol) = annihilation(Float64, B::Basis, i::Int, s)
annihilation(::Type{T}, B::Basis) where {T} = annihilation(T, B, i, :sparse)

"""
$(TYPEDSIGNATURES)
"""
creation(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = transpose(annihilation(T, B, i, s))
creation(B::Basis, i::Int, s::Symbol) = creation(Float64, B::Basis, i::Int, s)
creation(::Type{T}, B::Basis) where {T} = creation(T, B, i, :sparse)
=#