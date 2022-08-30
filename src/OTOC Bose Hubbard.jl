using .OTOC_example

L = 6
N = 6
basis1 = Szbasis(L,N)

basis2 = Szbasis(L,N-1)

basis3 = Szbasis(L,N-1)


H1 = sparse_hamiltonian(basis1, N)
      
H2 = sparse_hamiltonian(basis2, N-1)
H3 = sparse_hamiltonian(basis3, N-2)

t_vals = range(0,0.2,40)
vals = OTOC_lattice.(6,1,t_vals)

scatter(t_vals,vals)