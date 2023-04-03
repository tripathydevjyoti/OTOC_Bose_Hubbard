using LinearAlgebra, DifferentialEquations

# Define the differential equation with time-dependent coefficients for matrices
function f(du,u,p,t)
    a = time_dependent_coefficients(t, p) # call a user-defined function to get the coefficients
    du .= a*u*a # the differential equation
end

# Define a user-defined function to get the time-dependent coefficients
function time_dependent_coefficients(t, p)
    return [exp(1im* t*p[1]) -exp(-1im* t*p[1]); exp(1im* t*p[1]) exp(-1im* t*p[1])]
end

# Set the initial condition
u0 = [3.0 + 1im 0.0+0.5im; 0.0+0.5im 1.0-1im]

# Set the time span and time step
tspan = (0.0, 10.0)
dt = 0.1

# Set the parameter vector for the user-defined function
p = [1.0]

# Define the problem
prob = ODEProblem(f, u0, tspan, p)

# Solve the problem using the Tsit5 solver
sol = solve(prob, Tsit5(), dt=dt)

# Print the solution
println(sol)

using Plots

# Get the time and solution vectors from the solution object
t = sol.t
u = sol.u
u
print(t)
imag(u[2][1])

plot_u = []
for i in 1:length(t)
    push!(plot_u, real(u[i][1]))

end
plot_u
# Plot the variation of the 1st diagonal element of the solution matrix with respect to time
plot(t, plot_u, xlabel="Time", ylabel="u[1,1]", title="Solution Matrix")



