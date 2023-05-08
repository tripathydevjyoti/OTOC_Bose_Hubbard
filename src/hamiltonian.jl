export
    BoseHubbard

"""
$(TYPEDSIGNATURES)
"""
struct BoseHubbard{T <: Number}
    basis::Union{Basis, NBasis, Vector{NBasis},NplusBasis}
    lattice::LabelledGraph
    H::SparseMatrixCSC{T, Int64}
end

function BoseHubbard{T}(B, lattice::LabelledGraph, dummy::Int) where T <: Number
    H = spzeros(T, B.dim, B.dim)

    for edge ∈ edges(lattice)
        J = get_prop(lattice, edge, :J)
        H -= J .* annihilation(T, B, dst(edge))' * annihilation(T, B, src(edge))
    end
    H += H'

    for i ∈ vertices(lattice)
        U = get_prop(lattice, i, :U) ./ T(2)
        n = occupation(T, B, i)
        H += U .* (n * n - n)
    end
    BoseHubbard{T}(B, lattice, H)
end

function BoseHubbard{T}(B, lattice::LabelledGraph) where T <: Number
    I, J, V = Int[], Int[], T[]
    U = get_prop.(Ref(lattice), vertices(lattice), :U) ./ T(2)

    for (v, ket) ∈ enumerate(B.eig_vecs)
    # 1. interaction part:
        push!(I, v)
        push!(J, v)
        push!(V, sum(U .* ket .* (ket .- 1)))

    # 2. kinetic part:
        for edge ∈ edges(lattice)
            i, j = src(edge), dst(edge)
            if ket[j] > 0 && ket[i] != B.N
                w = get_index(B, create(destroy(ket, j), i))
                Jij = get_prop(lattice, edge, :J)
                push!(I, w, v)
                push!(J, v, w)
                val = -Jij * sqrt((ket[i] + 1) * ket[j])
                push!(V, val, conj(val))
            end
        end
    end
    BoseHubbard{T}(B, lattice, sparse(I, J, V, B.dim, B.dim))
end

BoseHubbard(B, lattice::LabelledGraph) = BoseHubbard{Float64}(B, lattice)

function BoseHubbard(N::IntOrVec, M::Int, J::T, U::T, bndr::Symbol) where T <: Number
    BoseHubbard{T}(NBasis(N, M), chain(M, J, U, bndr))
end

function BoseHubbard(B, J::T, U::T, graph) where T <: Number
    inst = Dict((src(e), dst(e)) => J for e ∈ edges(graph))
    push!(inst, ((i, i) => U for i ∈ 1:nv(graph))...)
    BoseHubbard{T}(B, lattice(T, inst))
end


function BoseHubbard(N::IntOrVec, M::Int, J::T, U::T, graph) where T <: Number
    BoseHubbard(NBasis(N, M), J, U, graph)
end


Base.eltype(ham::BoseHubbard{T}) where {T} = T
