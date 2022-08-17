export
    tag,
    Basis

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

       #basis =[[σ...] for σ ∈ Iterators.product(fill(collect(0:N), M)...)] |> vec

        tags = tag.(basis)
        order = sortperm(tags)
        new{Int, eltype(tags)}(N, M, tags[order], basis[order])
    end
end
