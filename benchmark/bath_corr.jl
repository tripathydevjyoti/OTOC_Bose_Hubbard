include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

using KrylovKit
using LinearAlgebra






function bath_bartek(dim::Dims, time1::Real ,time2::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time1)

    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    #times = range(0,1.2,50)
    state = State([one(T)], [fill(1, M)])
    bath(time1, time2, Ref(H), site1, site2, Ref(state))
end  

function bath_bartek2(dim::Dims, time1::Real, time2::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time1)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    #times = range(0,1.2,50)
    state = State([one(T)], [fill(1, M)])
    
    bath2(time1, time2, Ref(H), site1, site2, Ref(state))
end   



dim = (1, 1)
time1 =0.25
time2 = 0.1
num_points = 40

@time bath_corr11 = bath_bartek(dim,time1,time2,num_points,4,1)
