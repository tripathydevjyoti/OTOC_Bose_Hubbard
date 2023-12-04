include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard
using Test
using LightGraphs
using LabelledGraphs
using MetaGraphs
using Plots
using PyCall



function otoc_bartek(dim::Dims, time::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time)
    J, U = T(4), T(16)
    graph = twod_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0,2.0,40)
    state = State([one(T)], [fill(1, M)])
    OTOC.(times, Ref(H), site1, site2, Ref(state))
end    

dim = (1, 2)
time1 = 0
num_points = 40
@time otoc = otoc_bartek(dim, time1, num_points, 2,5)
abs.(otoc)
twothre = abs.(otoc)

Jtimes = range(0,2.0,40)
plot(Jtimes, [twothre,abs.(otoc)])
np = pyimport("numpy")
np.save("figfive_two_five.npy",abs.(otoc))
np.save("otoc_1hex_mi",otoc)
#print(otoc)

"""
nx = pyimport("networkx")
dim=(1,2)
bndr=false
hg = nx.generators.lattice.hexagonal_lattice_graph(
    dim..., periodic = bndr == :PBC ? true : false
)
map = Dict(v => i for (i, v) ∈ enumerate(hg.nodes))
print(map)
J=4
inst = Dict((map[v], map[w]) => J for (v, w) ∈ hg.edges)
print(inst)
push!(inst, ((map[v], map[v]) => U for v ∈ hg.nodes)...)
"""
"""
    g = cycle_graph(14)
    add_edge!(g,3,12)
    add_edge!(g,5,10)

    nodes1 = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2))
    map = Dict(v => i for (i, v) ∈ enumerate(nodes1))

    edges1 =  [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (0, 2)), ((0, 2), (1, 2)), ((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, 1), (2, 1)), ((1, 2), (1, 3)), ((1, 3), (2, 3)), ((2, 0), (2, 1)), ((2, 0), (3, 0)), ((2, 1), (2, 2)), ((2, 2), (2, 3)), ((2, 2), (3, 2)), ((3, 0), (3, 1)), ((3, 1), (3, 2))]
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ edges1)
    push!(inst, ((map[v], map[v]) => U for v ∈ nodes1)...)

    lattice(T, inst)


    function hex_graph(J::T, U::T)  where T <: Real
    
        g = cycle_graph(6)
        nodes1 = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2))
        map = Dict(v => i for (i, v) ∈ enumerate(nodes1))
    
        edges1 =  [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (0, 2)), ((0, 2), (1, 2)), ((1, 0), (1, 1)), ((1, 1), (1, 2))]
        inst = Dict((map[v], map[w]) => J for (v, w) ∈ edges1)
        
        push!(inst, ((map[v], map[v]) => U for v ∈ nodes1)...)
        lattice(T, inst)
    end  



"""

"""
function hexagonal_graph(dim::Dims, J::T, U::T, bndr::Symbol) where T <: Real
    @assert bndr ∈ (:OBC, :PBC)

    nx = pyimport("networkx")
    hg = nx.generators.lattice.hexagonal_lattice_graph(
        dim..., periodic = bndr == :PBC ? true : false
    )
    map = Dict(v => i for (i, v) ∈ enumerate(hg.nodes))
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ hg.edges)
    push!(inst, ((map[v], map[v]) => U for v ∈ hg.nodes)...)

    lattice(T, inst)
end
"""

"""
function hexagonal_graph(dim::Dims, J::T, U::T, bndr::Symbol) where T <: Real
    @assert bndr ∈ (:OBC, :PBC)

    nx = pyimport("networkx")
    hg = nx.generators.lattice.hexagonal_lattice_graph(
        dim..., periodic = bndr == :PBC ? true : false
    )
    map = Dict(v => i for (i, v) ∈ enumerate(hg.nodes))
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ hg.edges)
    push!(inst, ((map[v], map[v]) => U for v ∈ hg.nodes)...)

    lattice(T, inst)
end
"""

nx =pyimport("networkx")
g  =nx.Graph()
g.add_nodes_from([1,2,3,4,5,6,7,8,9,10,11,12,13])
g.nodes()
g.add_edges_from([(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(1,12),(4,13),(8,13)])