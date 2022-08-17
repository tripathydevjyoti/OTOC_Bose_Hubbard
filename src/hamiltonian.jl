export
    hamiltonian

"""
$(TYPEDSIGNATURES)
"""
function hamiltonian(::Type{T}, B::Basis, lattice::LabelledGraph) where T
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
                tg = tag(destroy_and_create(ket, i, j))
                w = searchsortedfirst(B.tags, tg)
                Jij = get_prop(lattice, edge, :J)
                push!(J, v, w)
                push!(I, w, v)
                val = -Jij * sqrt((ket[j] + 1) * ket[i])
                push!(V, val, conj(val))
            end
        end
    end
    sparse(I, J, V, n, n)
end
hamiltonian(B::Basis, lattice::LabelledGraph) = hamiltonian(Float64, B, lattice)

function hamiltonian(N::Int, M::Int, J::T, U::T, boundry::Symbol) where T <: Real
    hamiltonian(Basis(N, M), bose_bubbard_1D(M, J, U, Val(boundry)))
end
