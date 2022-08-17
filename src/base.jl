export
    Basis, operate,
    destroy_and_create,
    creation, annihilation, hamiltonian

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
    for (v, ket) ∈ eig_vecs, (w, bra) ∈ eig_vecs
        if ket[i] > 0
            if bra == operate(ket, i, :destroy)
                a[w, v] = T(ket[i]) |> sqrt
            end
        end
    end
    a
end

annihilation(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = annihilation(T, B, i, Val(s))
annihilation(B::Basis, i::Int, s::Symbol) = annihilation(Float64, B::Basis, i::Int, s)
annihilation(B::Basis, i::Int) = annihilation(B, i, :sparse)

"""
$(TYPEDSIGNATURES)
"""
creation(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = transpose(annihilation(T, B, i, s))
creation(B::Basis, i::Int, s::Symbol) = creation(Float64, B::Basis, i::Int, s)
creation(B::Basis, i::Int) = creation(B, i, :sparse)

"""
$(TYPEDSIGNATURES)
"""
function hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph, ::Val{:sparse}) where T
    n = length(B.eig_vecs)

    # 1. interaction part:
    I = [collect(1:n)]
    J = copy(I)
    U = get_prop.(Ref(lattice), I, Ref(:U)) ./ 2
    V = T[sum(U .* ket .* (ket .- 1)) for ket ∈ B.eig_vecs]

    # 2. kinetic part:
    for (v, ket) ∈ enumerate(B.eig_vecs)
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

    eig_vecs = enumerate(B.eig_vecs)
    for (v, ket) ∈ eig_vecs
        # 1. interaction part:
        U = get_prop.(Ref(lattice), vertices(lattice), Ref(:U)) ./ 2
        H[v, v] = sum(U .* ket .* (ket .- 1)) |> T

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

hamiltonian(::Type{T}, B::Basis, i::Int, s::Symbol) where {T} = hamiltonia(T, B, i, Val(s))
hamiltonian(B::Basis, i::Int, s::Symbol) = hamiltonia(Float64, B::Basis, i::Int, s)
hamiltonian(B::Basis, i::Int) = hamiltonia(B, i, :sparse)

function hamiltonian(N::Int, M::Int, J::T, U::T, boundry::Symbol) where T <: Real
    hamiltonian(Basis(N, M), bose_bubbard_1D(M, J, U, Val(boundry)))
end
