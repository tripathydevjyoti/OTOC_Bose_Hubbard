export
    Basis,
    occupation,
    hamiltonian

#tag(v::Vector{Int}) = sum(log.(100 .* collect(1:length(v)) .+ 3) .* v)
tag(v::Vector{Int}) = hash(v)

"""
$(TYPEDSIGNATURES)
"""
struct Basis{T, S}
    N::T
    M::T
    tags::Vector{S}
    eig_vecs::Vector{Vector{T}}

    function Basis(N, M)
        basis = map(
            A -> [sum(A .== i) for i ∈ 1:N],
            with_replacement_combinations(1:N, M)
        )
        tags = tag.(basis)
        order = sortperm(tags)
        new{Int, eltype(tags)}(N, M, tags[order], basis[order])
    end
end

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, B::Basis, ::Val{:dense}) where T <: Real
    n = length(B.eig_vecs)
    N = zeros(T, n, n)
    for (v, ket) ∈ enumerate(B.eig_vecs)
        N[v, v] = sum(ket)
    end
    N
end

function occupation(::Type{T}, B::Basis, ::Val{:sparse}) where T <: Real
    n = length(B.eig_vecs)
    I, V = Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        push!(I, v)
        push!(V, sum(ket))
    end
    sparse(I, I, V, n, n)
end

occupation(::Type{T}, B::Basis, s::Symbol) where {T} = occupation(T, B, Val(s))
occupation(B::Basis, s::Symbol) = occupation(Float64, B::Basis, s)
occupation(::Type{T}, B::Basis) where {T} = occupation(T, B, :sparse)

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, B::Basis, i::Int, ::Val{:dense}) where T
    n = length(B.eig_vecs)
    N = zeros(T, n, n)
    for (v, ket) ∈ enumerate(B.eig_vecs)
        N[v, v] = ket[i]
    end
    N
end

function occupation(::Type{T}, B::Basis, i::Int, ::Val{:sparse}) where T
    n = length(B.eig_vecs)
    I, V = Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        push!(I, v)
        push!(V, ket[i])
    end
    sparse(I, I, V, n, n)
end

occupation(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = occupation(T, B, i, Val(s))
occupation(B::Basis, i::Int, s::Symbol) = occupation(Float64, B::Basis, i, s)
occupation(B::Basis, i::Int) = occupation(T, B, i, :sparse)

"""
$(TYPEDSIGNATURES)
"""
function hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph, ::Val{:sparse}) where T
    n = length(B.eig_vecs)

    I, J, V = Int[], Int[], T[]
    U = get_prop.(Ref(lattice), I, Ref(:U)) ./ 2

    for (v, ket) ∈ enumerate(B.eig_vecs)
    # 1. interaction part:
        push!(I, v)
        push!(J, v)
        push!(V, sum( U .* ket .* (ket .- 1)))

    # 2. kinetic part:
        for edge ∈ edges(lattice)
            i, j = src(edge), dst(edge)
            if ket[i] > 0 && ket[j] != B.N
                tg = tag(destroy_and_create(ket, i, j))
                w = searchsortedfirst(B.tags, tg)
                J = get_prop(lattice, edge, :J)
                push!(J, v, w)
                push!(I, w, v)
                val = -J * sqrt((ket[j] + 1) * ket[i])
                push!(V, val, conj(val))
            end
        end
    end
    sparse(I, J, V, n, n)
end

function hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph, ::Val{:dense}) where T <: Real
    n = length(B.eig_vecs)

    H = zeros(T, n, n)
    U = get_prop.(Ref(lattice), vertices(lattice), Ref(:U)) ./ 2

    eig_vecs = enumerate(B.eig_vecs)
    for (v, ket) ∈ eig_vecs
        # 1. interaction part:
        H[v, v] = sum(U .* ket .* (ket .- 1))

        # 2. kinetic part:
        for (w, bra) ∈ eig_vecs
            for edge ∈ edges(lattice)
                i, j = src(edge), dst(edge)
                if ket[i] > 0 && ket[j] != B.N
                    if bra == destroy_and_create(ket, i, j)
                        J = get_prop(lattice, edge, :J)
                        a[w, v] = -J * sqrt((ket[j] + 1) * ket[i])
                        a[v, w] = a[w, v]
                    end
                end
            end
        end
    end
    H
end

function hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph, s::Symbol) where {T}
    hamiltonian(T, B, lattice, Val(s))
end

hamiltonian(B::Basis, lattice::LabelledGraph, s::Symbol) =
    hamiltonian(Float64, B, lattice, s)

function hamiltonian(N::Int, M::Int, J::T, U::T, boundry::Symbol) where T <: Real
    hamiltonian(Basis(N, M), bose_bubbard_1D(M, J, U, Val(boundry)), :sparse)
end
