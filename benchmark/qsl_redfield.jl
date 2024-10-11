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
using ExponentialUtilities


tcheck=1.1
T = eltype(tcheck)

J, U = T(4), T(8)
N=M=2
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



"""

N = 3
M = 3
J = 4.0 # hopping parameter (float values only)
U = 9.0
beta = 1.0

T = eltype(J)
np = pyimport("numpy")

t_stop = 40.0

num_points = 15
times = np.linspace(0, t_stop, num_points)

# Initialize Hamiltonians and eigenvectors
red_ham = [RBoseHubbard(N+1, M, J, U), RBoseHubbard(N, M, J, U)]
eigenvals, eigenvecs = eigen(Matrix(red_ham[2].H))
init_dm = thermal_dm(red_ham[2], beta)
"""
#subsys_size = red_ham[2].basis.dim
np = pyimport("numpy")
t_stop = 40.0
num_points = 15
times = np.linspace(0, t_stop, num_points)
# Initialize arrays for the nested integrals
nested_integral1 = zeros(ComplexF64, length(times))
nested_integral2 = zeros(ComplexF64, length(times))
nested_integral3 = zeros(ComplexF64, length(times))
nested_integral4 = zeros(ComplexF64, length(times))

n = N
# Perform the integration over t
for (i, t_val) in enumerate(times)
    t_range = np.linspace(0, t_val, 100)

    # Initialize arrays for integration at each t'
    integral_results1 = zeros(ComplexF64, length(t_range))
    integral_results2 = zeros(ComplexF64, length(t_range))
    integral_results3 = zeros(ComplexF64, length(t_range))
    integral_results4 = zeros(ComplexF64, length(t_range))

    for (k, t_prime) in enumerate(t_range)
        # Initialize s_range for the inner integral
        s_range = range(0, stop=t_prime, length=100)

        # Compute time evolution of density matrix
        τ = 1im * t_prime
        #time_evol_dm = exponential!(τ * Matrix(red_ham[2].H)) * init_dm * exponential!(-τ * Matrix(red_ham[2].H))
        #red_dm = partial_trace(time_evol_dm, subsys_size, N, M)
        #red_dm = thermal_dm(red_ham[2], beta)

        # Arrays for the function values over s
        f_vals1 = ComplexF64[]
        g_vals1 = ComplexF64[]
        f_vals2 = ComplexF64[]
        g_vals2 = ComplexF64[]

        for j in 1:length(s_range)
            """
            push!(f_vals1, two_time_corr(red_ham, eigenvecs, [t_prime, s_range[j]], red_dm)[1])
            push!(g_vals1, two_time_corr(red_ham, eigenvecs, [t_prime, s_range[j]], red_dm)[2])
            push!(f_vals2, two_time_corr(red_ham, eigenvecs, [t_prime, s_range[j]], red_dm)[1])
            push!(g_vals2, two_time_corr(red_ham, eigenvecs, [t_prime, s_range[j]], red_dm)[2])
            """
            push!(f_vals1, thermal_corr(beta, H,t_prime, s_range[j]))[1]
            push!(g_vals1, thermal_corr(beta, H,t_prime, s_range[j]))[2]
            push!(f_vals2, thermal_corr(beta, H,t_prime, s_range[j]))[1]
            push!(g_vals2, thermal_corr(beta, H,t_prime, s_range[j]))[2]
        end

        sum_int = zeros(ComplexF64, 4)
        for k = 0:n
            for j in 1:length(s_range)
                f_vals1[j] *= exp(1im * U * (k - 1) * s_range[j])
                g_vals1[j] *= exp(-1im * U * (k - 1) * s_range[j])
            end

            integral_sf1 = trapz(s_range, f_vals1)
            sum_int[1] += (k) * integral_sf1 * conj(integral_sf1)

            integral_sg1 = trapz(s_range, g_vals1)
            sum_int[3] += (k) * integral_sg1 * conj(integral_sg1)
        end

        for k = 0:n-1
            for j in 1:length(s_range)
                f_vals2[j] *= exp(1im * U * (k) * s_range[j])
                g_vals2[j] *= exp(-1im * U * (k) * s_range[j])
            end

            integral_sf2 = trapz(s_range, f_vals2)
            sum_int[2] += (k+1) * integral_sf2 * conj(integral_sf2)

            integral_sg2 = trapz(s_range, g_vals2)
            sum_int[4] += (k+1) * integral_sg2 * conj(integral_sg2)
        end

        # Store the nested integral result for each t'
        integral_results1[k] = real(sqrt(sum_int[1]))
        integral_results2[k] = real(sqrt(sum_int[2]))
        integral_results3[k] = real(sqrt(sum_int[3]))
        integral_results4[k] = real(sqrt(sum_int[4]))
    end

    # Compute the nested integral over t' using the trapezoidal rule
    nested_integral1[i] = trapz(t_range, integral_results1)
    nested_integral2[i] = trapz(t_range, integral_results2)
    nested_integral3[i] = trapz(t_range, integral_results3)
    nested_integral4[i] = trapz(t_range, integral_results4)

    println(i)
end

# Final plot for the nested integral
using Plots

plot(times, real((nested_integral2 + nested_integral4) * sqrt(n*(n+1)/2) + (nested_integral1 + nested_integral3) * sqrt(n*(n+3)/2)), xlabel="Time", ylabel="Nested Integral", title="Nested Two-Time Correlation Integral")
