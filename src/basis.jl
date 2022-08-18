export
    Basis,
    get_index

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

    function Basis(N::Int, M::Int; constraint::Symbol=:none)
        @assert constraint ∈ (:none, :conserved_particles)

        if constraint == :conserved_particles
            basis = map(
                A -> [sum(A .== i) for i ∈ 1:N],
                with_replacement_combinations(1:N, M)
            )
        elseif constraint == :none
            basis =[[σ...] for σ ∈ Iterators.product(fill(collect(0:N), M)...)] |> vec
        else
            throw(DomainError("Attempt to specify $constraint which is not allowed."))
        end

        tags = tag.(basis)
        order = sortperm(tags)
        new{Int, eltype(tags)}(N, M, tags[order], basis[order])
    end
end

get_index(B::Basis, ket::Vector{Int}) = searchsortedfirst(B.tags, tag(ket))
