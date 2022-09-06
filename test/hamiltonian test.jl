using Test
@testset "Checking if both hamiltonian functions match" begin
    basis = Szbasis(6,6)
    H1 = sparse_hamiltonian(basis,6)

    graph = cycle_graph(6)
    H2 = graph_hamiltonian(basis, graph)

@test Matrix(H1) ==  Matrix(H2)
end
