using OTOC_Bose_Hubbard
using LabelledGraphs
using LightGraphs
using MetaGraphs
using Test

my_tests = [
    "base.jl",
    #"lattice.jl",
]

for my_test ∈ my_tests
    include(my_test)
end
