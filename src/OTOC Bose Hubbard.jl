using .OTOC_example
using JeszenszkiBasis

L = 6
N = 6
"""
basis1 = Szbasis(L,N)

basis2 = Szbasis(L,N-1)

basis3 = Szbasis(L,N-2)


H1 = sparse_hamiltonian(basis1, N)
      
H2 = sparse_hamiltonian(basis2, N)
H3 = sparse_hamiltonian(basis3, N)
"""
t_vals = range(0,0.2,80)
vals = []

for t in t_vals
    push!( vals, OTOC_lattice([1,1,1,1,1,1], 6, 1, t, L, N))
    
end    

scatter(t_vals,vals)