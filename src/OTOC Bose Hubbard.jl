using .OTOC_example
using JeszenszkiBasis
using Plots
using LightGraphs





g = cycle_graph(10)
add_edge!(g,3,8)
#grid with two hexagons





t_vals = range(0,1.0,50)
vals = []
site1 = 10
site2 = 1
M = 10
J=4                            #no of particles 
for t in t_vals
    push!( vals, OTOC_graphlattice([1,1,1,1,1,1,1,1,1,1], site1, site2, t, g, M))
    print(t)
    
end


plot(J*t_vals ,vals)