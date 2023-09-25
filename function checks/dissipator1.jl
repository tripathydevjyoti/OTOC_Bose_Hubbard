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

function dissipator1(time1::Real, time2::Real, site1::Int, site2::Int, den_mat)
    T = eltype(time1)
    J, U = T(0), T(16)
    graph = system_graph(J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NplusBasis.([N-1, N, N+1], M), Ref(graph))
    basis_states = [[2,1],[3,0],[1,2],[0,3],[1,1],[2,0],[0,2],[1,0],[0,1]]
    diss_mat1 = zeros(9,9)
    for i in 1:length(basis_states)
        for j in 1:length(basis_states)
            state1= State([one(T)], [basis_states[i]])
            state2= State([one(T)], [basis_states[j]])

            diss_mat1[i][j] = diss_one_cre(time1, time2, H, site1, site2, state1, state2, den_mat) 
            diss_mat1[i][j] += -diss_two_cre(time1, time2, H, site1, site2, state1, state2, den_mat) 
        end
    end    

    return diss_mat1    
end

den_mat = zeros(9,9)
den_mat[1,1]=1.0+0.0im
dissipator1(0.3,0.1,1,2,den_mat)