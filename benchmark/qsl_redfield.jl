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
using Trapz

tcheck=1.1
T = eltype(tcheck)

J, U = T(4), T(16)
N=M=6
beta = 1.0
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
function part_func(beta, H)
    
    part_sum = 0
    for i in 1:length(eigenvals)
        part_sum = part_sum + exp(-beta*eigenvals[i])
        
    end
    part_sum
end        

partition_function = part_func(T(1), H)


function thermal_corr(beta, H, time1, time2)
    

    trsum1 = 0
    trsum2 = 0
   
    for i in 1:length(eigenvals)
        
        gibbs = exp(-beta*eigenvals[i])
        gamma1 = gibbs*bath(time1, time2, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        gamma2 = gibbs*bath2(time1, time2, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
       
        trsum1 = trsum1 + gamma1
        trsum2 = trsum2 + gamma2
        
    end
    
    corr_arr = [trsum1/ partition_function, trsum2/partition_function]
return corr_arr
end



# Define the upper limit for t
T = 0.1

# Define the range of integration for t
t_range = range(0, stop=T, length=100)

# Initialize array to store intermediate integration results
integral_results1 = zeros(length(t_range))
integral_results2 = zeros(length(t_range))
integral_results3 = zeros(length(t_range))
integral_results4 = zeros(length(t_range))


n=1
# Integrate over t
for (i, t_val) in enumerate(t_range)
    # Define the range of integration for s (from 0 to t_val)
    s_range = range(0, stop=t_val, length=100)
    
    f_vals1 = ComplexF64[]
    g_vals1 = ComplexF64[]
    f_vals2 = ComplexF64[]
    g_vals2 = ComplexF64[]

    for j in 1:length(s_range)
        push!(f_vals1, thermal_corr(beta, H, t_val, s_range[j])[1]) 
        push!(g_vals1, thermal_corr(beta, H, t_val, s_range[j])[2])
        push!(f_vals2, thermal_corr(beta, H, t_val, s_range[j])[1])
        push!(g_vals2, thermal_corr(beta, H, t_val, s_range[j])[2]) 
    end
    # Evaluate the function over the range of s
    sum_int = zeros(4)
    for k in 0:n 
        
      
        

        for j in 1:length(s_range)
            f_vals1[j] = f_vals1[j]*exp(1im*U*(k-1)*s_range[j])
            g_vals1[j] = g_vals1[j]*exp(-1im*U*(k-1)*s_range[j])

        end
   
        integral_sf1 = trapz(s_range, f_vals1)
        sum_int[1] = sum_int[1] + (k)*integral_sf1*conj(integral_sf1)
        
        integral_sg1 = trapz(s_range, g_vals1)
        sum_int[3] = sum_int[3] + (k)*integral_sg1*conj(integral_sg1)
    end    
  
    for k in 0:n-1
      

        for j in 1:length(s_range)
            f_vals2[j] = f_vals2[j]*exp(1im*U*(k)*s_range[j])
            g_vals2[j] = g_vals2[j]*exp(-1im*U*(k)*s_range[j])
        end

        integral_sf2 = trapz(s_range, f_vals2)
        sum_int[2] = sum_int[2] + (k+1)*integral_sf2*conj(integral_sf2)

        integral_sg2 = trapz(s_range, g_vals2)
        sum_int[4] = sum_int[4] + (k+1)*integral_sg2*conj(integral_sg2)
    
    end    
    # Store the result for this t_val
    integral_results1[i] = sum_int[1]
    integral_results2[i] = sum_int[2]
    integral_results3[i] = sum_int[3]
    integral_results4[i] = sum_int[4]
    print(i)
end

# Integrate the results over t using the trapezoidal rule
final_integral1 = trapz(t_range, integral_results1)
final_integral2 = trapz(t_range, integral_results2)
final_integral3 = trapz(t_range, integral_results3)
final_integral4 = trapz(t_range, integral_results4)

final_integral = (f_vals2 + g_vals2)*sqrt(n*(n+1)/2) + (f_vals1 + g_vals1)*sqrt(n*(n+3)/2)

println("Double integral result: ", final_integral)


