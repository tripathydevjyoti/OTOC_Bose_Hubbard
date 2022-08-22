export
    operate,
    create, destroy,
    creation, annihilation

"""
$(TYPEDSIGNATURES)
"""
function operate(ket::Vector{Int}, i::Int, c::Int)
    nket = copy(ket)
    nket[i] += c
    nket
end

create(ket::Vector{Int}, i::Int) = operate(ket, i, 1)
destroy(ket::Vector{Int}, i::Int) = ket[i] > 0 ? operate(ket, i, -1) : 0

function destroy(state::State{T}, i::Int) where T
    n = length(state.eig_vecs)
    vecs = Vector(undef, n)
    Threads.@threads for k ∈ 1:n
        vecs[k] = destroy(state.eig_vecs[k], i)
    end
    K = findall(iszero, vecs)
    State(state.coeff[K], state.eig_vecs[K])
end

"""
$(TYPEDSIGNATURES)
"""
function annihilation(::Type{T}, B::Basis, i::Int) where T <: Real
    n = length(B.eig_vecs)
    I, J, V = Int[], Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        if ket[i] > 0
            push!(J, v)
            push!(I, get_index(B, operate(ket, i, -1)))
            push!(V, T(ket[i]) |> sqrt)
        end
    end
    sparse(I, J, V, n, n)
end
annihilation(B::Basis, i::Int) = annihilation(Float64, B::Basis, i::Int)

"""
$(TYPEDSIGNATURES)
"""
creation(::Type{T}, B::Basis, i::Int) where {T} = transpose(annihilation(T, B, i))
