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



t_lower = 0.0
np = pyimport("numpy")
t_upper = np.linspace(0,1.0,12)

s_lower = 0.0
s_upper = 10.0
sprime_lower = 0.0
sprime_upper = 10.0


function integrand_f(t)
    T = eltype(t)

    J, U = T(4), T(16)
    N=M=5
    H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
    state = State([one(T)], [fill(1, M)])
    s_upper = t

    integrate_s, _= quadgk(s -> bath(t, s, H, 1, 1, state), s_lower, s_upper)
    integrate_s
end

function integrand_g(t)
    T = eltype(t)

    J, U = T(4), T(16)
    N=M=5
    H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC)
    state = State([one(T)], [fill(1, M)])
    sprime_upper = t

    integrate_sprime, _= quadgk(sprime -> bath2(t, sprime, H, 1, 1, state), sprime_lower, sprime_upper)
    integrate_sprime
end

@time result_f, _ = quadgk(integrand_f, t_lower, t_upper)
@time result_g, _ = quadgk(integrand_g, t_lower, t_upper)

result = Complex[]
for i in 1:length(t_upper)
    result_f, _ = quadgk(integrand_f, t_lower, t_upper[i])
    result_g, _ = quadgk(integrand_g, t_lower, t_upper[i])  
    push!(result, result_f*result_g)
end    

plot(t_upper, abs.(result))
np.save("qsl.npy",result)


