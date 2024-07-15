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



t_lower = 0.0
np = pyimport("numpy")
t_upper = np.linspace(0,1.0,12)
t_upper = 0.1

s_lower = 0.0
s_upper = 10.0
sprime_lower = 0.0
sprime_upper = 10.0


n=2
tcheck=1.1
T = eltype(tcheck)

J, U = T(4), T(16)
N=M=6
beta = 1.0
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
#eigenvalues, eigenstates = eigs(H[2].H, nev=1, which=:SR)
#lowest_eigenstate = eigenstates[:,1]
#state = State(lowest_eigenstate,H[2].basis)
#state = State([one(T)], [fill(1, M)])
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

function integrand_fdag(t)
    T = eltype(t)

    s_upper = t

    sum_int = 0
    for i in 0:n
        integrate_s, _= quadgk(s -> (bath(t, s, H, 1, 1, state)*exp(1im*U*(i-1)*s)), s_lower, s_upper)
        integrate_sprime, _= quadgk(s -> conj(bath(t, s, H, 1, 1, state)*exp(1im*U*(i-1)*s)), s_lower, s_upper)
        sum_int = sum_int + i*integrate_s*integrate_sprime
    end    
    sqrt(sum_int)
end

function integrand_f(t)
    T = eltype(t)

    
    s_upper = t

    sum_int = 0
    for i in 0:n-1
        integrate_s, _= quadgk(s -> (bath(t, s, H, 1, 1, state)*exp(1im*U*(i)*s)), s_lower, s_upper)
        integrate_sprime, _= quadgk(s -> conj(bath(t, s, H, 1, 1, state)*exp(1im*U*(i)*s)), s_lower, s_upper)
        sum_int = sum_int + (i+1)*integrate_s*integrate_sprime
    end    
    sqrt(sum_int)
end

function integrand_g(t)
    T = eltype(t)

    s_upper = t

    sum_int = 0
    for i in 0:n
        integrate_s, _= quadgk(s -> (bath2(t, s, H, 1, 1, state)*exp(-1im*U*(i-1)*s)), s_lower, s_upper)
        integrate_sprime, _= quadgk(s -> conj(bath2(t, s, H, 1, 1, state)*exp(-1im*U*(i-1)*s)), s_lower, s_upper)
        sum_int = sum_int + i*integrate_s*integrate_sprime
    end    
    sqrt(sum_int)
end


function integrand_gdag(t)
    T = eltype(t)


    s_upper = t

    sum_int = 0
    for i in 0:n-1
        integrate_s, _= quadgk(s -> (bath2(t, s, H, 1, 1, state)*exp(-1im*U*i*s)), s_lower, s_upper)
        integrate_sprime, _= quadgk(s -> conj(bath2(t, s, H, 1, 1, state)*exp(-1im*U*i*s)), s_lower, s_upper)
        sum_int = sum_int + (i+1)*integrate_s*integrate_sprime
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
function int_lind2(t, ns)
    s_upper = t
    integrate_s, _=quadgk(s -> (bath2(t, s, H, 1, 1, state))*exp(1im*U*ns*s), s_lower, s_upper)

    return integrate_s
end
"""


int_lind1(0.1,0)
t_upper = np.linspace(0,1.0,12)
result = []

#result_lind =[]
result_lind = int_lind1(0,2)
bound_sum = 0
for j in 0:n-1
    result_lind = int_lind1(0.1,j)
    bound_sum = bound_sum + (i+1)*(np.real((result_lind[1]) )+ np.real( (result_lind[2])))

end
print(bound_sum)



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
