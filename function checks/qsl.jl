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


function integrate_thermal_corr(beta, H, time_range)
    
    
    integral_trsum1 = trapz(time_range, trsum1_vals)
    print("int1done")
    integral_trsum2 = trapz(time_range, trsum2_vals)
    
    return integral_trsum1, integral_trsum2
end


J = 4.0
U = 16.0
N=M=6             #no of sites and bosons in the chain
beta = 1.0        #inverse temperature
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC) #BH hamiltonian 
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
eigenvals
partition_function = part_func(1.0, H)  #partition function for the BH hamiltonian


np = pyimport("numpy")
time_range = np.linspace(0,0.9,100)
trsum1_vals = Float64[]
trsum2_vals = Float64[]

    
@showprogress for t in time_range
        corr_arr = thermal_corr(beta, H, t,partition_function)
        push!(trsum1_vals, np.real(corr_arr[1]))
        push!(trsum2_vals, np.real(corr_arr[2]))
        print(t)
end
int_sum_lin2 =np.save("trsum2_vals.npy", trsum_vals)
integral_trsum1, integral_trsum2 = integrate_thermal_corr(beta, H, time_range)


n=1
bound_sum = 0
for j in 0:n-1
    
    bound_sum = bound_sum + (j+1)*(np.real((integral_trsum1) )+ np.real( (integral_trsum2)) )

end
print(bound_sum)


times = np.linspace(0,0.1,20)
plot(J*times, 4*J*J*bound_sum*times)
np.save("qsl_lind.npy", 4*J*J*bound_sum*times)
qsl = np.load("qsl_lind.npy")
plot(J*times, qsl)
"""

for i in 1:length(t_upper)
    result_f, _ = quadgk(integrand_f , t_lower, t_upper[i])
    result_fdag, _ = quadgk(integrand_fdag , t_lower, t_upper[i])
    result_g, _ = quadgk(integrand_g , t_lower, t_upper[i])
    result_gdag, _ = quadgk(integrand_gdag , t_lower, t_upper[i])  
    push!(result, (result_f + result_gdag)*sqrt(n*(n+1)/2) + (result_fdag+result_g)*sqrt(n*(n+3)/2))
end    

plot(t_upper, [4*J*J*(np.real(result)),  4*J*J*result_lind])
plot(t_upper, 4*16*result_lind)
plot(t_upper, 4*16*(np.real(result)) )

store =4*16*(np.real(result))
store_lin = 4*16*result_lind
np.save("qsl_check.npy",store)
np.save("qsl_lind.npy", store_lin)


function integrand_fdag(t)
    T = eltype(t)

    s_upper = t

    sum_int = 0
    for i in 0:n
        integrate_s, _= quadgk(s -> ( (thermal_corr(beta, H, t, s)[1]) *exp(1im*U*(i-1)*s)), s_lower, s_upper)
        #integrate_sprime, _= quadgk(s -> conj(bath(t, s, H, 1, 1, state)*exp(1im*U*(i-1)*s)), s_lower, s_upper)
        sum_int = sum_int + i*integrate_s*conj(integrate_s)
    end    
    sqrt(sum_int)
end

function integrand_f(t)
    T = eltype(t)

    
    s_upper = t

    sum_int = 0
    for i in 0:n-1
        integrate_s, _= quadgk(s -> ((thermal_corr(beta, H, t, s)[1])*exp(1im*U*(i)*s)), s_lower, s_upper)
        #integrate_sprime, _= quadgk(s -> conj(bath(t, s, H, 1, 1, state)*exp(1im*U*(i)*s)), s_lower, s_upper)
        sum_int = sum_int + (i+1)*integrate_s*conj(integrate_s)
    end    
    sqrt(sum_int)
end

function integrand_g(t)
    T = eltype(t)

    s_upper = t

    sum_int = 0
    for i in 0:n
        integrate_s, _= quadgk(s -> ((thermal_corr(beta, H, t, s)[2])*exp(-1im*U*(i-1)*s)), s_lower, s_upper)
        #integrate_sprime, _= quadgk(s -> conj(bath2(t, s, H, 1, 1, state)*exp(-1im*U*(i-1)*s)), s_lower, s_upper)
        sum_int = sum_int + i*integrate_s*conj(integrate_s)
    end    
    sqrt(sum_int)
end


function integrand_gdag(t)
    T = eltype(t)


    s_upper = t

    sum_int = 0
    for i in 0:n-1
        integrate_s, _= quadgk(s -> ((thermal_corr(beta, H, t, s)[2])*exp(-1im*U*i*s)), s_lower, s_upper)
        #integrate_sprime, _= quadgk(s -> conj(bath2(t, s, H, 1, 1, state)*exp(-1im*U*i*s)), s_lower, s_upper)
        sum_int = sum_int + (i+1)*integrate_s*conj(integrate_s)
    end    
    sqrt(sum_int)
end    
#@time result_f, _ = quadgk(integrand_f , t_lower, t_upper)
#@time result_g, _ = quadgk(integrand_g , t_lower, t_upper)

function int_lind1(t, ns)

   
    s_upper = 0.1
    integrate_s, _=quadgk(s -> (thermal_corr(1.0, H, s)[1])*exp(-1im*U*ns*s), s_lower, s_upper)
    integrate_sprime, _=quadgk(s -> (thermal_corr(1.0, H, s)[2])*exp(1im*U*ns*s), s_lower, s_upper)

    return [integrate_s, integrate_sprime]
end
"""


"""
function partitionft(
    beta::T, H::Vector{BoseHubbard{S}}; kwargs=()
) where{S, T <:Real}
    
    part_sum = 0
    for (i, comp_basis) in enumerate(H[2].basis.eig_vecs)
        state = State([1.0], [comp_basis])
        half_bra = expv( -beta, H[2], state)
        part_sum = part_sum + dot(dense(state,H[2].basis), half_bra)

    end
    
    return part_sum
end    



function expv(τ::Number, ham::BoseHubbard{T}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(
        ham.H, τ, dense(v, ham.basis), ishermitian=true, tol=1E-8
    )
    @assert info.converged == 1
    U_dket
end
"""