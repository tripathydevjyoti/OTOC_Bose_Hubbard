using .OTOC_example
using JeszenszkiBasis
using Plots




g = cycle_graph(10)
add_edge!(g,3,8)
#grid with two hexagons

t_vals = range(0,0.25,80)
vals = []
vals1 = []
for t in t_vals
    push!( vals, OTOC_graphlattice([1,1,1,1,1,1,1,1,1,1], 10, 1, t, g, 10))
    
end

for t in t_vals
    push!( vals1, OTOC_lattice([1,1,1,1,1,1], 6, 1, t, 6, 6))
end    

plot(t_vals, [vals,vals1])