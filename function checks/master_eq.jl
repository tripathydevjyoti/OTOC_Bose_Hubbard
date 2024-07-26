include("../src/ME_Bose_Hubbard.jl")
using .ME_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs
using DifferentialEquations
using Plots
using PyCall
using QuadGK
using KrylovKit
using LinearAlgebra
using Combinatorics


function rhs!(du, u, p, t)
    # Unpack the state variables
    ρs = reshape(u, (9, 9)) # Density matrix ρs is a 2x2 matrix
    diss_mat1 = p[1] # System operator A
    diss_mat2 = p[2]
    bath_corr1 = p[3] # Rate matrix Γ
    bath_corr2 = p[4]

    # Compute the time derivative of the density matrix
    du_mat = zeros(9, 9)
    for α in 1:2
        for β in 1:2
            integrand(tau) = (
                bath_corr1(t,tau,α,β)*( dissipator1(t, tau, α, β, ρs)
                )
                           + (
                 bath_corr2(t,tau,α,β)*( dissipator2(t, tau, α, β, ρs)
                            )
                ))
            integral, _ = quadgk(integrand, 0, 10)
            du_mat += integral
        end
    end
    du .= du_mat[:]
end

# Define the initial conditions and time span
u0 = [0.0 + 0im, 0.0 + 0im, 0.0 + 0im, 0.0 + 0im,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # Initial density matrix ρs is |0><0|
u0[19] = 1.0+0im
tspan = (0.0, 5.0)





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

function bath_corr2( time1::Real, time2::Real, site1::Int, site2::Int)
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
            sum = sum + bath2(time1, time2, H, i, j, state)
        end    
    end

    return sum
  
    #return bath2(time1,time2,H,site1,site2,state)  
end



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

function dissipator2(time1::Real, time2::Real, site1::Int, site2::Int, den_mat)
    T = eltype(time1)
    J, U = T(0), T(16)
    graph = system_graph(J::T, U::T, :OBC)
    
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NplusBasis.([N-1, N, N+1], M), Ref(graph))
    basis_states = [[2,1],[3,0],[1,2],[0,3],[1,1],[2,0],[0,2],[1,0],[0,1]]
    diss_mat2 = zeros(9,9)
    for i in 1:length(basis_states)
        for j in 1:length(basis_states)
            state1= State([one(T)], [basis_states[i]])
            state2= State([one(T)], [basis_states[j]])

            diss_mat2[i][j] = diss_one_des(time1, time2, H, site1, site2, state1, state2, den_mat) 
            diss_mat2[i][j] += -diss_two_des(time1, time2, H, site1, site2, state1, state2, den_mat) 
        end
    end    

    return diss_mat2    
end



p = [dissipator1, dissipator2, bath_corr1, bath_corr2]

# Define the ODE problem
prob = ODEProblem(rhs!, u0, tspan, p)

# Solve the ODE using the default ODE solver

sol = solve(prob)





"""
time1 =0.3
time2 = 0.1
num_points = 40

rho_t = [0.2+0im  0 0 0 0 0 0 0 0;
         0 0.2 0 0 0 0 0 0 0;
         0 0 0.1 0 0 0 0 0 0;
         0 0 0 0.1 0 0 0 0 0;
         0 0 0 0 0.1 0 0 0 0;
         0 0 0 0 0 0.1 0 0 0;
         0 0 0 0 0 0 0.1 0 0;
         0 0 0 0 0 0 0 0.05 0;
         0 0 0 0 0 0 0 0 0.05]

#note that this ρ(t) is not the correct density matrix. 
#as mentioned in the pdfs, we need to work with a hilbert space which can contain N-1,N,N+1 bosons or upto N+1 bosons.         




@time dissipator_11 = dissipator(time1,time2,1,1,rho_t)

@time bath_11 = bath_corr2(time1,time2,4,4)

"""



    

















