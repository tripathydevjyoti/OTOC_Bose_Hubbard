
module OTOC_Bose_Hubbard
    using LabelledGraphs
    using LightGraphs
    using MetaGraphs
    using CSV
    using KrylvKit
    using DocStringExtensions
    using LinearAlgebra, MKL

    include("lattice.jl")
    include("base.jl")
    include("OTOC.jl")

end # modul
