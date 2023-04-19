include("../src_compare/ME_Bose_Hubbard.jl")
using .ME_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

using KrylovKit
using LinearAlgebra








function bath_corr1( time1::Real ,time2::Real, site1::Int, site2::Int)
    T = eltype(time1)

    J, U = T(4), T(16)
    graph = hexagonal_graph( J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    
    state = State([one(T)], [fill(1, M)])
    bath(time1, time2, H, site1, site2, state)
end  

function bath_corr2( time1::Real, time2::Real, site1::Int, site2::Int)
    T = eltype(time1)
    J, U = T(4), T(16)
    graph = hexagonal_graph( J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    
    state = State([one(T)], [fill(1, M)])
    
    bath2(time1, time2, H, site1, site2, state)
end

function dissipator(time1::Real, time2::Real, site1::Int, site2::Int, den_mat)
    T = eltype(time1)
    J, U = T(0), T(16)
    graph = system_graph(J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1], M), Ref(graph))
  
    state1= State([one(T)], [fill(1, M)])
    state2= State([one(T)], [fill(1, M)])
    
    diss_one_cre(time1, time2, H, site1, site2, state1, state2,den_mat)
end
#note this is not the full dissiptor function



dim = (1, 1)
time1 =0.3
time2 = 0.1
num_points = 40

rho_t = [1 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0]

#note that this ρ(t) is not the correct density matrix. 
#as mentioned in the pdfs, we need to work with a hilbert space which can contain N-1,N,N+1 bosons or upto N+1 bosons.         

@time dissipator_11 = dissipator(time1,time2,1,1,rho_t)

@time bath_11 = bath_corr2(dim,time1,time2,1,1)













