include("../src_compare/ME_Bose_Hubbard.jl")
using .ME_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

using KrylovKit
using LinearAlgebra








function bath_bartek(dim::Dims, time1::Real ,time2::Real, site1::Int, site2::Int)
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
    bath(time1, time2, H, site1, site2, state)
end  

function bath_bartek2(dim::Dims, time1::Real, time2::Real, site1::Int, site2::Int)
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
    
    bath2(time1, time2, H, site1, site2, state)
end

function dissipator(time1::Real, time2::Real, site1::Int, site2::Int, den_mat)
    T = eltype(time1)
    J, U = T(0), T(16)
    graph = system_graph(J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    #times = range(0,1.2,50)
    state1= State([one(T)], [fill(1, M)])
    state2= State([one(T)], [fill(1, M)])
    
    diss_one_cre(time1, time2, H, site1, site2, state1, state2,den_mat)
end



dim = (1, 1)
time1 =0.3
time2 = 0.1
num_points = 40

rho_t = [1 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0]

@time dissipator_11 = dissipator(time1,time2,1,1,rho_t)

@time bath_11 = bath_bartek2(dim,time1,time2,1,1)













