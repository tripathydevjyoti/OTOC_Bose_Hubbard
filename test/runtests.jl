using OTOC_Bose_Hubbard
using LabelledGraphs
using LightGraphs
using MetaGraphs
using Test

my_tests = [
    "basis.jl",
    "lattice.jl",
    "conserved_operators.jl",
    "create_destroy.jl",
    "hamiltonian.jl",
]

for my_test ∈ my_tests
    include(my_test)
end
