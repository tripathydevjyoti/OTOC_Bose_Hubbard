export
    Basis,
    NBasis,
    get_index,
    dense_eigen_vec

abstract type AbstractBasis end

"""
$(TYPEDSIGNATURES)
"""
tag(v::Vector{Int}) = hash(v)

"""
$(TYPEDSIGNATURES)
"""
all_states(N::Int, M::Int) = vec(collect.(Iterators.product(fill(collect(0:N), M)...)))

"""
$(TYPEDSIGNATURES)
"""
function all_sub_states(N::Int, M::Int)
    basis = Vector{Int}[]
    η = [N, zeros(Int, M-1)...]
    push!(basis, η)

    for _ ∈ 1:M
        while η[1] > 0
            η = copy(η)
            η[1] -= 1
            η[2] += 1
            push!(basis, η)
        end
        i = findfirst(k->k != 0, basis[end])
        if i == M return break end
        η = copy(basis[end])
        η[i+1] += 1
        η[1] = η[i] - 1
        η[i] = 0
        push!(basis, η)
    end
    basis
end

struct Basis{T, S} <: AbstractBasis
    N::T
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}
    sub_basis_indices::Vector{T}

    function Basis(N::Int, M::Int)
        basis = all_states(N, M)
        tags = tag.(basis)
        order = sortperm(tags)
        B = basis[order]
        new{Int, eltype(tags)}(
            N, M, (N + 1) ^ M, tags[order], B, findall(v->sum(v) == N, B)
        )
    end
end

struct NBasis{T, S} <: AbstractBasis
    N::T
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function NBasis(N::Int, M::Int)
        basis = all_sub_states(N, M)
        tags = tag.(basis)
        order = sortperm(tags)
        new{Int, eltype(tags)}(
            N, M, length(basis), tags[order], basis[order]
        )
    end
end

"""
$(TYPEDSIGNATURES)
"""
get_index(B, ket::Vector{Int}) = searchsortedfirst(B.tags, tag(ket))

"""
$(TYPEDSIGNATURES)
"""
function dense_eigen_vec(::Type{T}, B, ket::Vector{Int}) where T <: Real
    dket = zeros(T, length(B.eig_vecs))
    dket[get_index(B, ket)] = one(T)
    dket
end
dense_eigen_vec(B, ket::Vector{Int}) = dense_eigen_vec(Float64, B, ket)
