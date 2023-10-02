export
    lattice,
    chain,
    hexagonal_graph,
    system_graph

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
function chain(M::Int, J::T, U::T, bndr::Symbol) where T <: Real
    @assert bndr ∈ (:OBC, :PBC)
    inst = Dict((i, i+1) => J for i ∈ 1:M-1)
    push!(inst, ((i, i) => U for i ∈ 1:M)...)
    if bndr == :PBC push!(inst, (M, 1) => J) end
    lattice(T, inst)
end



"""
$(TYPEDSIGNATURES)
"""





function hexagonal_graph( J::T, U::T, bndr::Symbol) where T <: Real
    @assert bndr ∈ (:OBC, :PBC)

    nx = pyimport("networkx")
    hg  = nx.Graph()
    hg.add_nodes_from([1,2,3,4,5,6,7,8])
    hg.add_edges_from([(1,2),(2,3),(4,5),(5,6),(6,7),(7,8)])

    map = Dict(v => i for (i, v) ∈ enumerate(hg.nodes))
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ hg.edges)
    push!(inst, ((map[v], map[v]) => U for v ∈ hg.nodes)...)

    lattice(T, inst)
end




"""
$(TYPEDSIGNATURES)
"""


function system_graph(J::T, U::T, bndr::Symbol) where T <: Real
    @assert bndr ∈ (:OBC, :PBC)

    nx = pyimport("networkx")
    hg  = nx.Graph()
    hg.add_nodes_from([1,2])
    hg.add_edges_from([(1,2)])

    map = Dict(v => i for (i, v) ∈ enumerate(hg.nodes))
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ hg.edges)
    push!(inst, ((map[v], map[v]) => U for v ∈ hg.nodes)...)

    lattice(T, inst)
end






