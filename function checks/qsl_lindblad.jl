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

np = pyimport("numpy")
ns=2
tcheck=1.1
T = eltype(tcheck)

J, U = T(4), T(16)
N=M=6
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
#eigenvalues, eigenstates = eigs(H[2].H, nev=1, which=:SR)
#lowest_eigenstate = eigenstates[:,1]
#state = State(lowest_eigenstate,H[2].basis)
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
function part_func(beta, H)
    
    part_sum = 0
    for i in 1:length(eigenvals)
        part_sum = part_sum + exp(-beta*eigenvals[i])
        
    end
    part_sum
end        




function thermal_corr(beta, H, time)
    

    #eigenvals = eigen_result.values
    #eigenvecs = eigen_result.vectors
    trsum1, trsum1_conj, trsum2, trsum2_conj = 0, 0, 0, 0
    for i in 1:length(eigenvals)
        gamma1 =bath(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis), beta )
        #gamma1_conj = conj(gamma1)
        gamma2 = bath2(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis), beta )
        #gamma2_conj = conj(gamma2)
        
        trsum1 = trsum1 + gamma1
        #trsum1_conj = trsum1_conj + gamma1_conj
        trsum2 = trsum2 + gamma2
        #trsum2_conj = trsum2_conj +gamma2_conj
    end
    partition_function = part_func(beta, H)
    corr_arr = [trsum1/ partition_function, trsum2/partition_function]
return corr_arr
end


arr = thermal_corr(T(10) , H, 0.1)

twopt1 = []
twopt1_conj = []
twopt2 = []
twopt2_conj = []

beta = T(10)
times = np.linspace(0, 2.0 , 15)
for i in 1:length(times)
    arr = thermal_corr(beta, H, times[i])
    print(arr)
    push!(twopt1, arr[1])
    #push!(twopt1_conj, arr[2])
    push!(twopt2, arr[2])
    #push!(twopt2_conj, arr[4])
end    

plot(times, [real(twopt1),  real(twopt2)])