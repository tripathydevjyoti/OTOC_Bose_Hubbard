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


function thermal_corr(                                                   #function to compute thermal correlators
    beta::T, H::Vector{BoseHubbard{S}}, time::T, partition::T;  kwargs=()
    ) where{S, T <:Real}
    
    trsum1 = 0
    trsum2 = 0
    
    eigenvals, eigenvecs = eigen(Matrix(H[2].H))
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



np = pyimport("numpy")
timearg = np.linspace(0, 0.01, 5) #array for times
coupling = np.linspace(0, 6, 6)   #array for range of coupling U/J

F1real =zeros( length(coupling), length(timearg) ) # Matrix to store real part of Gamma1
F1imag =zeros( length(coupling), length(timearg) ) # Matrix to store imag part of Gamma1
F2real =zeros( length(coupling), length(timearg) ) # Matrix to store real part of Gamma2
F2imag =zeros( length(coupling), length(timearg) ) # Matrix to store imag part of Gamma2 



@showprogress for i in 1:length(coupling)                     #loop to store values for heatmap
    J = 4.0
    U = coupling[i]*J
    N=M=6 
    beta = 1.0
    H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
    
    partition_function = part_func(beta, H)
    
    for j in 1:length(timearg)
        
        
        twoptarr =  thermal_corr(beta, H, timearg[j], partition_function)
        
        F1real[i,j] = np.real(twoptarr[1])
        F1imag[i,j] = np.imag(twoptarr[1])
        F2real[i,j] = np.real(twoptarr[2])
        F2imag[i,j] = np.imag(twoptarr[2])
        print(timearg[j])
    end
    
end


heatmaps = [
    (F1real, "F1real", "ℜ(Γ1)"),
    (F1imag, "F1imag", "ℑ(Γ1)"),
    (F2real, "F2real", "ℜ(Γ2)"),
    (F2imag, "F2imag", "ℑ(Γ2)")
]

for (matrix, filename_suffix, title_label) in heatmaps
    heatmap(timearg, coupling, matrix, xlabel="Time", ylabel="Coupling", title=title_label, color=:viridis)
    savefig(@sprintf("colormap_%s.pdf", title_label))
    
end
