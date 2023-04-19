export
    occupation

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, B) where T <: Real
    I = Vector{Int}(undef, B.dim)
    V = Vector{T}(undef, B.dim)
    for v ∈ 1:B.dim
        I[v] = v
        V[v] = sum(B.eig_vecs[v])
    end
    sparse(I, I, V, B.dim, B.dim)
end

"""
$(TYPEDSIGNATURES)
"""
function occupation(::Type{T}, B, i::Int) where T
    I = Vector{Int}(undef, B.dim)
    V = Vector{T}(undef, B.dim)
    for v ∈ 1:B.dim
        I[v] = v
        V[v] = B.eig_vecs[v][i]
    end
    sparse(I, I, V, B.dim, B.dim)
end
occupation(B) = occupation(Float64, B)
occupation(B, i::Int) = occupation(Float64, B, i)
