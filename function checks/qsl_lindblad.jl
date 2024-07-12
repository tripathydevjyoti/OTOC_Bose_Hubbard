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
state = State([one(T)], [fill(1, M)])
bath(0.01, 0.01, H, 1, 1, state)
function integrand1(x,n)
    return exp(-1im*U*n*x)*bath(x, x, H, 1, 1, state)
end

function integrand2(x,n)
    return exp(1im*U*n*x)*bath2(x, x, H, 1, 1, state)
end




sum =0
for i in 0:(ns-1)
    result1, error1 = quadgk(x -> integrand1(x, i), 0, 3)
    result2, error2 = quadgk(x -> integrand2(x, i), 0, 3)
    sum = sum + (i+1)*(np.real(result1)+np.real(result2))
    print(result1)
end

time = np.linspace(0, 0.1,12)
plot(time, 8*16*sum*time)
    
