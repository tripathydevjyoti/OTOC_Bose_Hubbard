export
    Basis,
    NBasis,
    State,
    get_index,
    RBasis,
    tensor_basis,
    limit_tensor_basis,
    dense

abstract type AbstractBasis end

const IntOrVec = Union{Int, Vector{Int}}

"""
$(TYPEDSIGNATURES)
"""
tag(v::Vector{T}) where T = hash(v)
tag(v) = hash(v)

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

struct RBasis{T, S} <: AbstractBasis
    N::Union{T, Vector{T}}
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function RBasis(states::Vector, N::IntOrVec, M::Int)
        # Generate tags for each state but don't reorder
        tags = tag.(states)
        
        # Keep the original order of states and tags
        new{Int, eltype(tags)}(N, M, length(states), tags, states)
    end

    RBasis(N::IntOrVec, M::Int) = RBasis(vcat(vcat.([NBasis(i, M-1).eig_vecs for i in N:-1:0])...), N, M-1)
end

"""
$(TYPEDSIGNATURES)
"""

struct tensor_basis{T, S} <: AbstractBasis
    N::Union{T, Vector{T}}
    M::T
    dim::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function tensor_basis(states, N::IntOrVec, M::Int)
        tags = tag.(states)

        new{Int, eltype(tags)}(N, M, length(states), tags, states)
    end
    
    function tensor_basis(N::Int, M::Int)
        sys_basis = RBasis(N,2).eig_vecs
        bath_basis = RBasis(N,M).eig_vecs
        products  = collect.(Iterators.product(sys_basis,bath_basis))
        #states = vec(transpose([vcat(p...) for p in products]))
        states = vec(permutedims([vcat(p...) for p in products]))
       

        return tensor_basis(states, N, M)

    end

    
end 

"""
$(TYPEDSIGNATURES)
"""
function limit_tensor_basis(N::Int, M::Int)
    og_states = tensor_basis(N,M).eig_vecs
    filtered_states = filter(state -> sum(state) == N, og_states)

# Return a new tensor_basis instance with the filtered states
    return tensor_basis(filtered_states, N, M)
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
end

Base.eltype(state::State{T}) where {T} = T

"""
$(TYPEDSIGNATURES)
"""
@inline get_index(B::Union{NBasis, Basis}, ket::Vector{Int})  = searchsortedfirst(B.tags, tag(ket))
@inline get_index(B::RBasis, ket::Vector{Int64}) = findfirst(x ->x ==tag(ket), B.tags)




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


"""

"""