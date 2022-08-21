export
    Basis,
    NBasis,
    State,
    get_index,
    dense

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
all_sub_states(N::Int, M::Int) = collect(multiexponents(M, N))

"""
$(TYPEDSIGNATURES)
"""
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

"""
$(TYPEDSIGNATURES)
"""
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
struct State{T}
    coeff::Vector{T}
    eig_vecs::Vector{Vector{Int}}

    State(coeff, vecs) = new{eltype(coeff)}(coeff, vecs)
end

Base.eltype(state::State{T}) where T = eltype(eltype(state.coeff))

"""
$(TYPEDSIGNATURES)
"""
get_index(B::T, ket::Vector{Int}) where T <: AbstractBasis = searchsortedfirst(B.tags, tag(ket))

"""
$(TYPEDSIGNATURES)
"""
function dense(::Type{T}, ket::Vector{Int}, B::S) where {T <: Real, S <: AbstractBasis}
    dket = zeros(T, B.dim)
    dket[get_index(B, ket)] = one(T)
    dket
end
dense(ket::Vector{Int}, B::T) where T <: AbstractBasis = dense(Float64, ket, B)

function dense(state::State, B)
     sum(state.coeff .* dense.(Ref(eltype(state)), state.eig_vecs, Ref(B)))
end
