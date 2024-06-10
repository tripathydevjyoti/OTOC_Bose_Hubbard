export
    Basis,
    NBasis,
    NplusBasis,
    State,
    get_index,
    dense

abstract type AbstractBasis end

const IntOrVec = Union{Int, Vector{Int}}

"""
$(TYPEDSIGNATURES)
"""
tag(v::Vector{T}) where T = hash(v)

"""
$(TYPEDSIGNATURES)
"""
all_states(N::Int, M::Int) = vec(collect.(Iterators.product(fill(collect(0:N), M)...)))


"""
$(TYPEDSIGNATURES)
"""
all_sub_states(N::Int, M::Int) = collect(multiexponents(M, N))
all_sub_states(N::Vector{Int}, M::Int) = vcat(all_sub_states.(N, M)...)



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
    N::Union{T, Vector{T}}
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function NBasis(states::Vector, N::IntOrVec, M::Int)
        tags = tag.(states)
        order = sortperm(tags)
        new{Int, eltype(tags)}(N, M, length(states), tags[order], states[order])
    end
    NBasis(N::IntOrVec, M::Int) = NBasis(all_sub_states(N, M), N, M)
end

"""
$(TYPEDSIGNATURES)
"""

struct NplusBasis{T, S} <: AbstractBasis
    N::Union{T, Vector{T}}
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function NplusBasis(states::Vector, N::IntOrVec, M::Int)
        tags = tag.(states)
        order = sortperm(tags)
        new{Int, eltype(tags)}(N, M, length(states), tags[order], states[order])
    end
    NplusBasis(N::IntOrVec, M::Int) = NplusBasis(all_sub_states([N-1,N,N+1], M), N, M)
end




"""
$(TYPEDSIGNATURES)
"""
struct State{T}
    coeff::Vector{T}
    eig_vecs::Vector{Vector{Int}}

    function State(coeff::Vector{T}, vecs::Vector) where T <: Number
        new{eltype(coeff)}(coeff, vecs)
    end

    function State(ket::Vector{T}, B::S) where {T <: Number, S <: AbstractBasis}
        K = findall(!iszero, ket)
        State(ket[K], B.eig_vecs[K])
    end

    function State(ket::Vector{T}, B::S, C::Vector{Any}) where {T <: Number, S <: AbstractBasis}
        K = findall(!iszero, ket)
        common_vecs = intersect(B.eig_vecs[K],C)
        new_indices =[]
        for i in 1:length(common_vecs)
            for j in 1:length(C)

                if common_vecs[i] == C[j]
                    append!(new_indices,j)
                end 
            end       
        end


        State(ket[K], C[new_indices])
    end
end


Base.eltype(state::State{T}) where {T} = T

"""
$(TYPEDSIGNATURES)
"""
@inline get_index(B::T, ket::Vector{Int}) where T <: AbstractBasis = searchsortedfirst(B.tags, tag(ket))

"""
$(TYPEDSIGNATURES)
"""
function dense(::Type{T}, eket::Vector{Int}, B::S) where {T <: Number, S <: AbstractBasis}
    dket = zeros(T, B.dim)
    dket[get_index(B, eket)] = one(T)
    dket
end
dense(ket::Vector{Int}, B::T) where T <: AbstractBasis = dense(Float64, ket, B)

function dense(state::State, B::T) where T <: AbstractBasis
    dket = zeros(eltype(state.coeff), B.dim)
    Threads.@threads for i ∈ 1:length(state.eig_vecs)
        dket[get_index(B, state.eig_vecs[i])] = state.coeff[i]
    end
    dket
end

