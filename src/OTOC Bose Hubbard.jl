using .OTOC_example
using JeszenszkiBasis

L = 6
N = 6

t_vals = range(0,1.0,80)
vals = []

for t in t_vals
    push!( vals, OTOC_lattice([1,1,1,1,1,1], 6, 1, t, L, N))
    
end    

scatter(t_vals,vals)