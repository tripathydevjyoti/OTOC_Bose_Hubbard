export
    Basis,
    creation, annihilation,
    bose_hubbard_hamiltonian

"""
$(TYPEDSIGNATURES)
"""
struct Basis{T <: Unsigned}
    N::T
    M::T
    tags::Vector{<:Unsigned}
    eig_vecs::Vector{Vector{T}}

    function Basis(N, M)
        basis = map(
            A -> [sum(A .== i) for i ∈ 1:n],
            with_replacement_combinations(1:n, k)
        )
        tags = hash.(basis)
        order = sortperm(tags)
        new{UInt}(N, M, tags[order], basis[order])
    end
end

"""
$(TYPEDSIGNATURES)
"""
function destroy!(ket::Vector{Int}, i::Int)
    if ket[i] == 0 return ket end
    nket = copy(ket)
    nket[i] -= 1
    nket
end

"""
$(TYPEDSIGNATURES)
"""
function create(ket::Vector{Int}, i::int)
    if ket[i] == sum(ket) return ket end
    nket = copy(ket)
    nket[i] += 1
    nket
end

"""
$(TYPEDSIGNATURES)
"""
destroy_and_create(ket::Vector{Int}, i::Int, j::Int) = create(destroy(ket, i), j)

"""
$(TYPEDSIGNATURES)
"""
function creation(B::Basis, i::Int)
    V, I, J = Int[], Int[], Float64[]
    for v ∈ B.tags
        ket = create(B.eig_vecs[v], i)
        w = searchsortedfirst(B.tags, hash(ket))
        push!(I, v)
        push!(J, w)
        push(V, sqrt((ket[i])))
    end
    sparse(V, I, J)
end

"""
$(TYPEDSIGNATURES)
"""
annihilation(B::Basis, i::Int) = transpose(creation(B, i))

"""
$(TYPEDSIGNATURES)
"""
function bose_hubbard_hamiltonian(B::Basis, lattice::LabelledGraph)
    V, I, J = Int[], Int[], Float64[]
    # 1. interaction part:

    # 2. kinetic part:
    for v ∈ B.tags
        for edge ∈ edges(lattice)
            i, j = src(edge), dst(edge)
            ket = destroy_and_create(ket, i, j)
            w = searchsortedfirst(B.tags, hash(ket))
            J = get_prop(lattice, edge, :J)
            push!(I, v, w)
            push!(J, w, v)
            val = J * sqrt((ket[i] + 1) * ket[j])
            push!(V, val, conj(val))
        end
    end
    sparse(V, I, J)
end
