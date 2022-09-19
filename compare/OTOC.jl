using .OTOC_Bose_Hubbard
using .OTOC_example
using Test
using LightGraphs
using LabelledGraphs

using Plots


function otoc_dev(site1::Int, site2::Int)
    g = cycle_graph(6)
    #add_edge!(g,3,8)
    
    t_vals = range(0,0.25,30)
    vals = []
    
    M = nv(g)
    
                        #no of particles 
    for t in t_vals
        push!( vals, OTOC_graphlattice([1,1,1,1,1,1], site1, site2, t, g, M))
    end
    vals
end
    

function otoc_bartek(dim::Dims, time::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    
    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0, 0.25, 30)
    state = State([one(T)], [fill(1, M)])
    
    J .* times, OTOC.(times, Ref(H), site1, site2, Ref(state))
end

dim = (1, 1)
time = 0.25
num_points = 40
Jtimes, otoc = otoc_bartek(dim, time, num_points, 3,1)
Jtimes, otoc1 = otoc_bartek(dim, time, num_points, 5,1)


#p1 =plot(Jtimes, [otoc_dev(6,1),otoc_dev(2,1)],label = ["OTOC(6,1)" "OTOC(2,1)"])
#p2 =plot(Jtimes, [otoc_dev(5,1),otoc_dev(3,1)],label = ["OTOC(5,1)" "OTOC(3,1)"])
p3 = plot(Jtimes, [abs.(otoc),abs.(otoc1)], label = ["OTOC(3,1)" "OTOC(5,1)"])
#p = plot(Jtimes, [abs.(otoc),otoc_dev(4,1)])
xlabel!("Jt")
ylabel!("OTOC(3,1) vs OTOC(5,1)")
savefig(p3, "./compare/otoc_bartek(3,1)vs(5,1).pdf")

"""
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
