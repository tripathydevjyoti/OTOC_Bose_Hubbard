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

partition_function = part_func(T(1), H)


function thermal_corr(beta, H, time)
    

    trsum1 = 0
    trsum2 = 0
   
    for i in 1:length(eigenvals)
        
        gibbs = exp(-beta*eigenvals[i])
        gamma1 = gibbs*bath(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        gamma2 = gibbs*bath2(time, time, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
       
        trsum1 = trsum1 + gamma1
        trsum2 = trsum2 + gamma2
        
    end
    
    corr_arr = [trsum1/ partition_function, trsum2/partition_function]
return corr_arr
end

@time bath(10.0, 10.0, H, 1, 1, State(eigenvecs[:, 1], H[2].basis)) 
@time arr = thermal_corr(1.0 , H, 10.0)
thermal_corr(1.0, H, 0.5)[1]
  


twopt1 = []
twopt1_conj = []
twopt2 = []
twopt2_conj = []

beta = T(1)
times = np.linspace(0, 1.0 , 30)
for i in 1:length(times)
    arr = thermal_corr(beta, H, times[i])
    print(arr)
    push!(twopt1, arr[1])
    #push!(twopt1_conj, arr[2])
    push!(twopt2, arr[2])
    #push!(twopt2_conj, arr[4])
end    

plot(times,  [ real(twopt1), real(twopt2), imag(twopt1), imag(twopt2)])
store =[]
for i in 1:length(twopt1)
    if i>15
        push!(store, real(twopt1[i]))
    end
end        
plot(times, store) 
