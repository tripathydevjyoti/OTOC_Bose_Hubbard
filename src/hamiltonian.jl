export
    BoseHubbard

"""
$(TYPEDSIGNATURES)
"""
struct BoseHubbard{T <: Real}
    basis
    lattice::LabelledGraph
    H::SparseMatrixCSC{T, Int64}
end

function BoseHubbard{T}(B, lattice::LabelledGraph) where T <: Real
    I, J, V = Int[], Int[], T[]
    U = get_prop.(Ref(lattice), vertices(lattice), Ref(:U)) ./ T(2)

    for (v, ket) ∈ enumerate(B.eig_vecs)
    # 1. interaction part:
        push!(I, v)
        push!(J, v)
        push!(V, sum(U .* ket .* (ket .- 1)))

    # 2. kinetic part:
        for edge ∈ edges(lattice)
            i, j = src(edge), dst(edge)
            if ket[i] > 0 && ket[j] != B.N
                w = get_index(B, create(destroy(ket, i), j))
                Jij = get_prop(lattice, edge, :J)
                push!(I, v, w)
                push!(J, w, v)
                val = -Jij * sqrt((ket[j] + 1) * ket[i])
                push!(V, val, conj(val))
            end
        end
    end
    BoseHubbard{T}(B, lattice, sparse(I, J, V, B.dim, B.dim))
end
BoseHubbard(B, lattice::LabelledGraph) = BoseHubbard{Float64}(B, lattice)

function BoseHubbard(N::Int, M::Int, J::T, U::T, bndr::Symbol) where T <: Real
    BoseHubbard{T}(NBasis(N, M), chain(M, J, U, bndr))
end

function BoseHubbard(B, J::T, U::T, graph) where T <: Real
    inst = Dict((src(e), dst(e)) => J for e ∈ edges(graph))
    push!(inst, ((i, i) => U for i ∈ 1:nv(graph))...)
    BoseHubbard{T}(B, lattice(T, inst))
end

function BoseHubbard(N::Int, M::Int, J::T, U::T, graph) where T <: Real
    BoseHubbard(NBasis(N, M), J, U, graph)
end
