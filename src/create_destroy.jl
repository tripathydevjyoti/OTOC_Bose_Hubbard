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
destroy(ket::Vector{Int}, i::Int) = operate(ket, i, -1)

function destroy(state::State{T}, i::Int) where T
    n = length(state.eig_vecs)
    vecs = Vector(undef, n)
    Threads.@threads for k ∈ 1:n
        ket = state.eig_vecs[k]
        vecs[k] = ket[i] > 0 ? operate(ket, i, -1) : 0
    end
    K = findall(!iszero, vecs)
    kets = vecs[K]
    State(state.coeff[K] .* sqrt.(getindex.(kets, i) .+ 1), kets)
end

function create(state::State{T}, i::Int) where T
    n = length(state.eig_vecs)
    vecs = Vector(undef, n)
    Threads.@threads for k ∈ 1:n
        ket = state.eig_vecs[k]
        vecs[k] =  operate(ket, i, +1)
    end
    K = findall(!iszero, vecs)
    kets = vecs[K]
    State(state.coeff[K] .* sqrt.(getindex.(kets, i) ), kets)
end

"""
$(TYPEDSIGNATURES)
"""
function annihilation(::Type{T}, B::S, i::Int) where {T <: Number, S <: AbstractBasis}
    n = length(B.eig_vecs)
    I, J, V = Int[], Int[], T[]
    for (v, ket) ∈ enumerate(B.eig_vecs)
        if ket[i] > 0
            push!(J, v)
            push!(I, get_index(B, operate(ket, i, -1)))
            push!(V, sqrt(T(ket[i])))
        end
    end
    sparse(I, J, V, n, n)
end
annihilation(B::Basis, i::Int) = annihilation(Float64, B::Basis, i::Int)

"""
$(TYPEDSIGNATURES)
"""
creation(::Type{T}, B::Basis, i::Int) where {T} = transpose(annihilation(T, B, i))
