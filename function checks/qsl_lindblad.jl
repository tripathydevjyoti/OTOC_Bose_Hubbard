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



J = 4.0           #hopping paramter
U = 8.0          #on-site potential
N=M=6             #no of sites and bosons in the chain
beta = 1.0        #inverse temperature
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC) #BH hamiltonian 
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
eigenvals
partition_function = part_func(1.0, H)  #partition function for the BH hamiltonian


times = np.linspace(0, 40.0, 50)
twopt1 =[]
twopt2 =[]



@showprogress for (_, t) in enumerate(times)
    arr = thermal_corr(beta, H, t, partition_function)
    push!(twopt1, arr[1])
    push!(twopt2, arr[2])
end    

#plot(times,  [ real(twopt1), real(twopt2), imag(twopt1), imag(twopt2)], label =["ℜ[Γ1]" "ℜ[Γ2]" "ℑ[Γ1]" "ℑ[Γ2]"], linewidth =2)
#plot(times, real(twopt1))
#xlabel!("t")
#title!("2-pt correlators for N=L=6 with β=1.0, U=4, J=4")
#savefig("6_6_BH_beta1_U4_J4_t5.pdf")
np.save("6_6_BH_beta1_U8_J4_t40(80)_2pt1.npy",twopt1)
np.save("6_6_BH_beta1_U8_J4_t40(80)_2pt2.npy",twopt2)
 

"""
figsee = np.load("6_6_BH_beta1_U8_J4_t10_2pt1.npy") 
plot(range(0,40,30),real(figsee))
fig1 = np.load("6_6_BH_beta1_U0_J4_t1_2pt1.npy")
fig2 = np.load("6_6_BH_beta1_U0_J4_t1_2pt2.npy")
fig3 = np.load("6_6_BH_beta1_U0_J4_t5_2pt1.npy")
fig4 = np.load("6_6_BH_beta1_U0_J4_t5_2pt2.npy")
fig5 = np.load("6_6_BH_beta1_U4_J4_t10_2pt1.npy")
fig6 = np.load("6_6_BH_beta1_U4_J4_t10_2pt2.npy")


times1 = np.linspace(0,1.0,30)
times2 = np.linspace(0,5.0,30)
times3 = np.linspace(0,10.0,30)
p1 = plot(times1,  [ real(fig1), real(fig2), imag(fig1), imag(fig2)], label =["ℜ[Γ1]" "ℜ[Γ2]" "ℑ[Γ1]" "ℑ[Γ2]"], linewidth =2)
p2 = plot(times2,  [ real(fig3), real(fig4), imag(fig3), imag(fig4)], label =["ℜ[Γ1]" "ℜ[Γ2]" "ℑ[Γ1]" "ℑ[Γ2]"], linewidth =2)
p3 = plot(times3,  [ real(fig5), real(fig6), imag(fig5), imag(fig6)], label =["ℜ[Γ1]" "ℜ[Γ2]" "ℑ[Γ1]" "ℑ[Γ2]"], linewidth =2)


plot(p1,p2, layout=(2,1), title="U/J=0, J=4, β=1.0", xlabel="t")
savefig("U0J4.pdf")
"""


