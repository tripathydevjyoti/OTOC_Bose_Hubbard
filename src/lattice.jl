export
    lattice,
    bose_bubbard_1D

    const Instance = Union{String, Dict}

"""
$(TYPEDSIGNATURES)
"""
function unique_nodes(bh_tuples)
    sort(collect(Set(Iterators.flatten((i, j) for (i, j, _) ∈ bh_tuples))))
end

"""
$(TYPEDSIGNATURES)
"""
function lattice(::Type{T}, instance::Instance) where T
    if instance isa String
        bhi = CSV.File(instance, types = [Int, Int, T], header=0, comment = "#")
    else
        bhi = [(i, j, J) for ((i, j), J) ∈ instance]
    end
    bhg = LabelledGraph{MetaGraph}(unique_nodes(bhi))

    set_prop!.(Ref(bhg), vertices(bhg), :U, T(0))

    for (i, j, v) ∈ bhi
        if i == j
            set_prop!(bhg, i, :U, T(v))
        else
            add_edge!(bhg, i, j) || throw(ArgumentError("Duplicate Egde ($i, $j)"))
            set_prop!(bhg, i, j, :J, T(v))
        end
    end
    bhg
end
lattice(instance::Instance) = lattice(Float64, instance)

"""
$(TYPEDSIGNATURES)
"""
function bose_bubbard_1D(M::Int, J::T, U::T, ::Val{:OBC}) where T <: Real
    inst = Dict{NTuple{2, Int}, T}((i, i+1) => J for i ∈ 1:M-1)
    push!(inst, ((i, i) => U for i ∈ 1:M)...)
    lattice(T, inst)
end
