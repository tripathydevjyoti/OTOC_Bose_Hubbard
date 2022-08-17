using OTOC_Bose_Hubbard
using LabelledGraphs
using LightGraphs
using MetaGraphs
using Test

my_tests = [
    "lattice.jl",
    "model.jl",
]

for my_test ∈ my_tests
    include(my_test)
end