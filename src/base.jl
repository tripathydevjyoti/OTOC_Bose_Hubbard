export
    Basis, operate,
    destroy_and_create,
    creation, annihilation,
    bose_hubbard_hamiltonian

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
function annihilation(::Type{T}, B::Basis, i::Int) where T
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
annihilation(B::Basis, i::Int) = annihilation(Float64, B::Basis, i::Int)

"""
$(TYPEDSIGNATURES)
"""
creation(B::Basis, i::Int) = transpose(annihilation(B, i))

"""
$(TYPEDSIGNATURES)
"""
function bose_hubbard_hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph) where T
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
                val = J * sqrt((ket[j] + 1) * ket[i])
                push!(V, val, conj(val))
            end
        end
    end
    sparse(I, J, V, n, n)
end
bose_hubbard_hamiltonian(B::Basis, lattice::LabelledGraph) =
bose_hubbard_hamiltonian(Float64, B, lattice)
