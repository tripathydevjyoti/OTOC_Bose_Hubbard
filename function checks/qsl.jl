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
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
eigenvalues, eigenstates = eigs(H[2].H, nev=1, which=:SR)
lowest_eigenstate = eigenstates[:,1]
#state = State(lowest_eigenstate,H[2].basis)
state = State([one(T)], [fill(1, M)])
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
    s_upper = t
    integrate_s, _=quadgk(s -> (bath(t, s, H, 1, 1, state))*exp(-1im*U*ns*s), s_lower, s_upper)

    return integrate_s
end

function int_lind2(t, ns)
    s_upper = t
    integrate_s, _=quadgk(s -> (bath2(t, s, H, 1, 1, state))*exp(1im*U*ns*s), s_lower, s_upper)

    return integrate_s
end


t_upper = np.linspace(0,0.1,12)
result = []

result_lind =[]


for i in 1:length(t_upper)
    sum =0
    
    for j in 0:n-1
        result1, _ = quadgk(t -> int_lind1(t, j), t_lower, t_upper[i] )
        result2, _ = quadgk(t -> int_lind2(t, j), t_lower, t_upper[i] )
        sum = sum + (i+1)*(np.real(result1)+ np.real(result2))
    end
    push!(result_lind, sum)
end



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
