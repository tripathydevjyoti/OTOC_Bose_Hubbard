using .OTOC_Bose_Hubbard
using .OTOC_example
using Test
using LightGraphs
using Plots


@testset "Compare both branches" begin
    function otoc_dev()
        g = cycle_graph(6)
        #add_edge!(g,3,8)
        
        t_vals = range(0,0.25,40)
        vals = []
        site1 = 6
        site2 = 1
        M = nv(g)
                            #no of particles 
        for t in t_vals
            push!( vals, OTOC_graphlattice([1,1,1,1,1,1], site1, site2, t, g, M))
        end
        vals
    end    

    function otoc_bartek(dim::Dims, time::Real, num_points::Int)
        T = eltype(time)
        J, U = T(4), T(16)
        graph = hexagonal_graph(dim, J::T, U::T, :OBC)
        M = nv(graph)
        N = Int(M / 1)
        H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

        #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
        times = range(0, 0.25, 40)
        state = State([one(T)], [fill(1, M)])
        J .* times, OTOC.(times, Ref(H), 2, 1, Ref(state))
    end

    dim = (1, 1)
    time = 0.25
    num_points = 40
    Jtimes, otoc = otoc_bartek(dim, time, num_points)

    plot(Jtimes, [abs.(otoc), otoc_dev()])

    @test isapprox(otoc_dev(), abs.(otoc) ) 


    
end






"""
function otoc_dev()
    g = cycle_graph(6)
    #add_edge!(g,3,8)
    #grid with two hexagons
    t_vals = range(0,0.7,40)
    vals = []
    
    site2 = 1
    M = nv(g)
    site1 = M
                        #no of particles 
    for t in t_vals
        push!( vals, OTOC_graphlattice([1,1,1,1,1,1], site1, site2, t, g, M))
    end
    vals
end    

function otoc_bartek(dim::Dims, time::Real, num_points::Int)
    T = eltype(time)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0, 0.25, 40)
    state = State([one(T)], [fill(1, M)])
    J .* times, OTOC.(times, Ref(H), 2, 1, Ref(state))
end

dim = (1, 1)
time = 0.25
num_points = 40
Jtimes, otoc = otoc_bartek(dim, time, num_points)

plot(Jtimes, [abs.(otoc), otoc_dev()])
"""


