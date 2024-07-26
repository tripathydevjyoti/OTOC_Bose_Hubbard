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
using Arpack
using ProgressMeter
using Printf



np = pyimport("numpy")
function thermal_corr(
    beta::T, H::Vector{BoseHubbard{S}}, time::T, partition::T;  kwargs=()
    ) where{S, T <:Real}
    
    trsum1 = 0
    trsum2 = 0
    
   
    for (i, energy) in enumerate(eigenvals)
        
        gibbs = exp(-beta*energy)
        
        gamma1 = gibbs*bath(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        trsum1 = trsum1 + gamma1 

        gamma2 = gibbs*bath2(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        trsum2 = trsum2 + gamma2
       
    end
    
    corr_arr = [trsum1/ partition, trsum2/partition ] 
return corr_arr
end



J = 4.0           #hopping paramter (float values only)
U = 8.0           #on-site potential (float values only)
N=M=6             #no of sites and bosons in the chain
beta = 1.0        #inverse temperature
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC) #BH hamiltonian 
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
eigenvals
partition_function = part_func(beta, H)  #partition function for the BH hamiltonian

t_stop = 5.0
num_points = 50
times = np.linspace(0, t_stop, num_points) #range of time values 


twopt1 =[] #array to store values for Γ1
twopt2 =[] #array to store values for Γ2



@showprogress for (_, t) in enumerate(times)
    arr = thermal_corr(beta, H, t, partition_function)
    push!(twopt1, arr[1])
    push!(twopt2, arr[2])
end    


filename1 = @sprintf("N_%d_L_%d_BH_beta_%.1f_U_%.1f_J_%.1f_t_%.1f_Gamma1.npy", N, M, beta, U, J, t_stop)
filename2 = @sprintf("N_%d_L_%d_BH_beta_%.1f_U_%.1f_J_%.1f_t_%.1f_Gamma2.npy", N, M, beta, U, J, t_stop)
np.save(filename1,twopt1)
np.save(filename2,twopt2)
 
 






