using OTOC_bose_bubbard
using Test

my_tests = [
    "base.jl",
    "lattice.jl",
]

for my_test ∈ my_tests
    include(my_test)
end
