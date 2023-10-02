include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard
using Test
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

dim=(1,1)
time1 =0.1
T = eltype(time1)
J, U = T(4), T(16)
graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
M = nv(graph)
N = Int(M / 1)
H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))
H[1].H

i=1
state = State([one(T)], [fill(1, M)])
length(H[1].basis.eig_vecs)