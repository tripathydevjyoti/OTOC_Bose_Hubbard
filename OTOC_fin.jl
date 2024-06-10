include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard
using Test
using LightGraphs
using LabelledGraphs

using Plots


function hex_graph(dim, J::T, U::T)
    g = cycle_graph(14)
    add_edge!(g,3,12)
    add_edge!(g,5,10)

    nodes1 = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2))
    map = Dict(v => i for (i, v) ∈ enumerate(nodes1))

    edges1 =  [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (0, 2)), ((0, 2), (1, 2)), ((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, 1), (2, 1)), ((1, 2), (1, 3)), ((1, 3), (2, 3)), ((2, 0), (2, 1)), ((2, 0), (3, 0)), ((2, 1), (2, 2)), ((2, 2), (2, 3)), ((2, 2), (3, 2)), ((3, 0), (3, 1)), ((3, 1), (3, 2))]
    inst = Dict((map[v], map[w]) => J for (v, w) ∈ edges1)
    push!(inst, ((map[v], map[v]) => U for v ∈ nodes1)...)

    lattice(T, inst)
end    


function otoc_bartek(dim::Dims, time::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)

    M = nv(graph)
    N = Int(M / 1)

    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0, 1.25, 50)
    state = State([one(T)], [fill(1, M)])
    J .* times, OTOC.(times, Ref(H), site1, site2, Ref(state))
end

dim = (1, 1)
time = 0.25
num_points = 40
Jtimes, otoc = otoc_bartek(dim, time, num_points, 2,1)
np = pyimport("numpy")
np.save("otoc_1hex_2_1",otoc)
#print(otoc)