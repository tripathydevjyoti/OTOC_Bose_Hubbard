
module OTOC_Bose_Hubbard
    using LabelledGraphs
    using LightGraphs
    using MetaGraphs
    using CSV
    using KrylovKit
    using SparseArrays
    using Combinatorics
    using DocStringExtensions
    using LinearAlgebra, MKL

    include("lattice.jl")
    include("base.jl")
    #include("OTOC.jl")

end # modul
