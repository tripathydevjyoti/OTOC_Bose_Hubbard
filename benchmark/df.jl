using LinearAlgebra
using DifferentialEquations
using QuadGK
include("../src/ME_Bose_Hubbard.jl")
using .ME_Bose_Hubbard
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

using KrylovKit
using LinearAlgebra
using Combinatorics



# Define the time derivative function
function rhs!(du, u, p, t)
    # Unpack the state variables
    ρs = reshape(u, (9, 9)) # Density matrix ρs is a 2x2 matrix
    A = p[1] # System operator A
    Adag = p[2]
    Γ = p[3] # Rate matrix Γ

    # Compute the time derivative of the density matrix
    du_mat = zeros(9, 9)
    for α in 1:2
        for β in 1:2
            integrand(tau) = (
                Γ[α, β](t-tau) * (
                    A[β](t-tau) * ρs * Adag[α](t) -
                    Adag[α](t) * A[β](t-tau) * ρs
                )
            )
            integral, _ = quadgk(integrand, 0, 10)
            du_mat += integral
        end
    end
    du .= du_mat[:]
end

# Define the initial conditions and time span
u0 = [1.0 + 0im, 0.0 + 0im, 0.0 + 0im, 0.0 + 0im,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # Initial density matrix ρs is |0><0|
tspan = (0.0, 5.0)

# Define the system operator A and rate matrix Γ
#A1(t) = [1 0; 0 exp(-im*t)] # A_1(t) operator
#
time1 = 0.3
T = eltype(time1)
J, U = T(0), T(16)
graph = system_graph(J::T, U::T, :OBC)

M = nv(graph)
N = Int(M / 1)
H = BoseHubbard.(NplusBasis.([N-1, N, N+1], M), Ref(graph))





a1dag = [ 0 0 sqrt(2) 0 0 0 0 0 0;
       0 0 0 1 0 0 0 0 0;
       0 0 0 0 0 sqrt(3) 0 0 0;
       0 0 0 0 0 0 sqrt(2) 0 0;
       0 0 0 0 0 0 0 1 0;
       0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0]

a2dag = [ 0 0 0 1 0 0 0 0 0;
          0 0 0 0 sqrt(2) 0 0 0 0;
          0 0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 0 sqrt(2) 0;
          0 0 0 0 0 0 0 0 sqrt(3);
          0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0]

a1 = [ 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       sqrt(2) 0 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       0 0 sqrt(3) 0 0 0 0 0 0;
       0 0 0 sqrt(2) 0 0 0 0 0;
       0 0 0 0 1 0 0 0 0;
       0 0 0 0 0 0 0 0 0] 

a2 = [ 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       1 0 0 0 0 0 0 0 0
       0 sqrt(2) 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0 0;
       0 0 0 sqrt(2) 0 0 0 0 0;
       0 0 0 0 sqrt(3) 0 0 0 0]       
       

A1(t) = exp(1im*t*Matrix(H[2].H)) * a1 * exp(-1im*t*Matrix(H[2].H))
A2(t) = exp(1im*t*Matrix(H[2].H)) * a2 * exp(-1im*t*Matrix(H[2].H))
A1dag(t) = exp(1im*t*Matrix(H[2].H)) * a1dag * exp(-1im*t*Matrix(H[2].H))
A2dag(t) = exp(1im*t*Matrix(H[2].H)) * a2dag * exp(-1im*t*Matrix(H[2].H))






#A2(t) = [0 exp(im*t); 0 0] # A_2(t) operator
A = [A1, A2] # Vector of system operators
Adag = [A1dag, A2dag]
Γ = [t -> (2*sin(t)+1)  t -> (2*t+1); t -> (2*t+1) t -> (2*t+1)] # Rate matrix

# Define the parameter vector
p = [A, Adag, Γ]

# Define the ODE problem
prob = ODEProblem(rhs!, u0, tspan, p)

# Solve the ODE using the default ODE solver

sol = solve(prob)

t = sol.t
ρ12 = [sol[i][2] for i in 1:length(t)]
using Plots
plot(t,imag(ρ12))






