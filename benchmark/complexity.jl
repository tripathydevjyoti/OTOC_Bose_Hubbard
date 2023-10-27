include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard
using Test
using LightGraphs
using LabelledGraphs
using LinearAlgebra
using Plots
using PyCall


time1 =0.1
T = eltype(time1)
J, U = T(4), T(10)
dim =(1,1)
graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
M = nv(graph)
N = Int(M / 1)
H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))
H[1].H

site =3
state = State([one(T)], [fill(1, M)])


dim_ham = length(H[1].basis.eig_vecs)
op_zero = zeros((dim_ham,dim_ham))
H[1].basis.eig_vecs

for i in 1:dim_ham
    op_zero[i,i] = float(H[1].basis.eig_vecs[i][site])
end    


function super_op(H,op)
    return H*op - op*H

end  



function norm(op)
    sqrt(tr(adjoint(op)*op))
end 

function trace_norm(op1,op2)
    tr(adjoint(op1)*op2)
end


op_zero = op_zero/norm(op_zero)

op_one = super_op(H[1].H,op_zero)
b1 = norm(op_one)
op_one = op_one/b1

op_basis = [op_one]
coeff_basis = [b1]

n_iter = 100

for i in 2:n_iter

    if i == 2
        new_op = super_op(H[1].H,op_basis[1]) - coeff_basis[1]*op_zero
        new_coeff = norm(new_op)
        push!(op_basis, new_op/new_coeff)
        push!(coeff_basis, new_coeff)

    else

        new_op = super_op(H[1].H,op_basis[i-1]) - coeff_basis[i-1]*op_basis[i-2]
        new_coeff = norm(new_op)
        push!(op_basis, new_op/new_coeff)
        push!(coeff_basis, new_coeff)
    end    

end    

coeff_basis
op_basis
#y_axis = [coeff_basis[1:10]]

x_axis = collect(1:n_iter)
sf_phase = coeff_basis
plot(x_axis,[coeff_basis,sf_phase], linewidth=3,  label=["U/J =4" "U/J=0.25"])
xlabel!("n")
ylabel!("Lanczos Coefficients")
savefig("/Users/dev/Desktop/coeff_BH.pdf")


Matrix(H[1].H)
function time_evol_op(time, op, hamil)
    exp(1im*time*hamil)*op*exp(-1im*time*hamil)
end

times = range(0.01, stop=5, length=20)
time_axis = collect(times)
complexity = []

for i in 1:length(time_axis)
    op_t = time_evol_op(time_axis[i], op_zero, Matrix(H[1].H))
    sum = 0
    phi_zero = trace_norm(op_t,op_zero)
    sum =sum + 0*abs(phi_zero)*abs(phi_zero)
    
    for j in 1:n_iter
        phi_j = trace_norm(op_t,op_basis[j])
        sum = sum + j*abs(phi_j)*abs(phi_j)
    end
    push!(complexity,sum)
end
sf_complexity = complexity
plot(time_axis,[mott_complexity,complexity],linewidth = 2.5,label =["U/J=4" "U/J=0.25"])
vline!([0.25], label="OTOC Decay", line=(color="black", linestyle=:dash))
xlabel!("t")
ylabel!("Krylov Complexity")
savefig("/Users/dev/Desktop/complexity_BH.pdf")



using LinearAlgebra
using Plots


# Compute eigenvalues
eigenvalues = eigvals(Matrix(H[1].H))

# Sort the eigenvalues in ascending order
sorted_eigenvalues = sort(eigenvalues)

# Calculate level spacings
level_spacings = diff(sorted_eigenvalues)
using Statistics
avg = mean(level_spacings)
# Plot the level spacing distribution
# Plot the level spacing distribution
histogram(level_spacings, bins=20, label="Level Spacing Distribution", legend=true)

# Customize the plot as needed
title!("Level Spacing Distribution for Quantum Chaos")
xlabel!("Level Spacing")
ylabel!("Frequency")

# Show the plot
plot()



    


