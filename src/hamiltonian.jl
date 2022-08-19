export
    BoseHubbard

"""
$(TYPEDSIGNATURES)
"""
struct BoseHubbard{T <: Real}
    basis::Basis
    lattice::LabelledGraph
    H::SparseMatrixCSC{T, Int64}
end

function BoseHubbard{T}(B::Basis, lattice::LabelledGraph) where T <: Real
    n = length(B.eig_vecs)

    I, J, V = Int[], Int[], T[]
    U = get_prop.(Ref(lattice), vertices(lattice), Ref(:U)) ./ 2

    for (v, ket) ∈ enumerate(B.eig_vecs)
    # 1. interaction part:
        push!(I, v)
        push!(J, v)
        push!(V, sum(U .* ket .* (ket .- 1)))

    # 2. kinetic part:
        for edge ∈ edges(lattice)
            i, j = src(edge), dst(edge)
            if ket[i] > 0 && ket[j] != B.N
                w = get_index(B, destroy_and_create(ket, i, j))
                Jij = get_prop(lattice, edge, :J)
                push!(I, v, w)
                push!(J, w, v)
                val = -Jij * sqrt((ket[j] + 1) * ket[i])
                push!(V, val, conj(val))
            end
        end
    end
    BoseHubbard{T}(B, lattice, sparse(I, J, V, n, n))
end
BoseHubbard(B::Basis, lattice::LabelledGraph) = BoseHubbard{Float64}(B, lattice)

function BoseHubbard(N::Int, M::Int, J::T, U::T, bndr::Symbol) where T <: Real
    BoseHubbard{T}(
        Basis(M, N; constraint=:conserved_particles),
        chain(M, J, U, Val(bndr))
    )
end
