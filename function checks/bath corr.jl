include("../src/ME_Bose_Hubbard.jl")
using .ME_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs
using QuadGK
using Plots
using PyCall
using DifferentialEquations
using KrylovKit
using LinearAlgebra
using Combinatorics


function bath_corr1( time1::Real ,time2::Real, site1::Int, site2::Int)
    T = eltype(time1)

    J, U = T(4), T(16)
    graph = hexagonal_graph( J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    
    state = State([one(T)], [fill(1, M)])


    if site1 == 1
        near_neib_site1 = [1,4]
    else
        near_neib_site1 = [3,4,8]    
    end
    
    if site2 == 1
        near_neib_site2 = [1,4]
    else
        near_neib_site2 = [3,4,8]    
    end
    
    sum = 0
    for i in near_neib_site1
        for j in near_neib_site2
            sum = sum + bath(time1, time2, H, i, j, state)
        end    
    end

    return sum
    
    
    #return bath(time1,time2,H,site1,site2,state)
end  


bath_corr1(0.3, 0.1, 1, 2)