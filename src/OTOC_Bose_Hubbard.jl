
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
    using DifferentialEquations
    using PyCall
    using CUDA, CUDA.CUSPARSE


    CUDA.allowscalar(false)

    include("basis.jl")
    include("lattice.jl")
    include("create_destroy.jl")
    include("conserved_operators.jl")
    include("hamiltonian.jl")
    include("OTOC.jl")

    # experimental:
    include("./experimental/OTOC_ODE.jl")
    include("./experimental/OTOC_CUDA.jl")


end # modul
