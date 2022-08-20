export
    occupation

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, eig_vecs::Vector{Vector{Int}}) where T <: Real
    dim = length(eig_vecs)
    I = Vector{Int}(undef, dim)
    V = Vector{T}(undef, dim)
    for v ∈ 1:dim
        I[v] = v
        V[v] = sum(eig_vecs[v])
    end
    sparse(I, I, V, dim, dim)
end

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, eig_vecs::Vector{Vector{Int}}, i::Int) where T
    dim = length(eig_vecs)
    I = Vector{Int}(undef, dim)
    V = Vector{T}(undef, dim)
    for v ∈ 1:dim
        I[v] = v
        V[v] = eig_vecs[v][i]
    end
    sparse(I, I, V, dim, dim)
end
occupation(eig_vecs::Vector{Vector{Int}}) = occupation(Float64, eig_vecs)
occupation(eig_vecs::Vector{Vector{Int}}, i::Int) = occupation(Float64, eig_vecs, i)
