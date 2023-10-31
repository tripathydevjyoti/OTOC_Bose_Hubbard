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
J, U = T(4), T(16)
dim =(1,1)
graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
M = nv(graph)
N = Int(M / 1)
H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))
#H = BoseHubbard(6,6 ,J, U, :OBC)
bh_hamil = H[1].H

site =3



dim_ham = length(H[1].basis.eig_vecs)
op_zero = zeros((dim_ham,dim_ham))
H[1].basis.eig_vecs

for i in 1:dim_ham
    op_zero[i,i] = float(H[1].basis.eig_vecs[i][site])
end    

op_zero

function super_op(H,op)
    return H*op - op*H

end  


function trace_norm(op1,op2)
    tr(adjoint(op1)*op2)
end

function norm(op)
    sqrt(trace_norm(op,op))
end 




op_zero = op_zero/norm(op_zero)

op_one = super_op(Matrix(bh_hamil),op_zero)
b1 = norm(op_one)
op_one = op_one/b1

op_basis = [op_one]
coeff_basis = [b1]

n_iter = 461
   

for i in 2:n_iter

    if i == 2
        new_op = super_op(bh_hamil,op_basis[1]) - coeff_basis[1]*op_zero
        new_coeff = norm(new_op)
        push!(op_basis, new_op/new_coeff)
        push!(coeff_basis, new_coeff)

    else

        new_op = super_op(bh_hamil,op_basis[i-1]) - coeff_basis[i-1]*op_basis[i-2]
        new_coeff = norm(new_op)
        push!(op_basis, new_op/new_coeff)
        push!(coeff_basis, new_coeff)
    end    

end    

coeff_basis
op_basis
#y_axis = [coeff_basis[1:10]]

x_axis = collect(1:n_iter)
coup4= coeff_basis
coup3 = coeff_basis
coup2 =coeff_basis
coup1 = coeff_basis
couphalf = coeff_basis
plot(x_axis,[couphalf,coup1,coup2,coup3,coup4,coeff_basis], linewidth=3,  label=["U/J =0.5" "U/J =1" "U/J =2" "U/J =3" "U/J =4"])
plot(x_axis,coeff_basis)
xlabel!("n")
ylabel!("Lanczos Coefficients")
savefig("/Users/dev/Desktop/coeff_BH.pdf")


Matrix(bh_hamil)
function time_evol_op(time, op, hamil)
    exp(1im*time*hamil)*op*exp(-1im*time*hamil)
end

times = range(0.0, stop=15, length=20)
time_axis = collect(times)
complexity = []

for i in 1:length(time_axis)
    op_t = time_evol_op(time_axis[i], op_zero,Matrix(bh_hamil))
    sum = 0
    phi_zero = trace_norm(op_zero,op_t)
    sum =sum + 0*abs(phi_zero)*abs(phi_zero)
    
    for j in 1:n_iter
        phi_j = trace_norm(op_basis[j],op_t)
        sum = sum + j*abs(phi_j)*abs(phi_j)
    end
    push!(complexity,sum)
end

plot(time_axis,[complexity],linewidth = 2.5,label =["U/J=4"])
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
num = length(level_spacings)
r_alpha = []
for i in 1:num
    if i!= num
        ra = level_spacings[i+1]/level_spacings[i]
        push!(r_alpha,min(ra,(1/ra)))
    end
end        
r_alpha       
histogram(r_alpha, bins=100 )
mean(r_alpha)

sort(level_spacings)
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



    


