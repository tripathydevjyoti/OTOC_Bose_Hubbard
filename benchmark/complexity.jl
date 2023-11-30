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


#H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))
H = BoseHubbard(5, 3 ,J, U, :OBC)
#bh_hamil = H.H
bh_hamil = H.H

site =2



dim_ham = length(H.basis.eig_vecs)
op_zero = zeros((dim_ham,dim_ham))
H.basis.eig_vecs


for i in 1:dim_ham
    op_zero[i,i] = (H.basis.eig_vecs[i][site])
end    

op_zero
adjoint(op_zero)==op_zero

function super_op(H,op)
    return H*op - op*H

end  


function trace_norm(op1,op2)
    tr(adjoint(op1)*op2)/dim_ham
end

function norm(op)
    sqrt(trace_norm(op,op))
end 




op_zero = op_zero/norm(op_zero)
norm(op_zero)
op_one = super_op(Matrix(bh_hamil),op_zero)
b1 = norm(op_one)
op_one = op_one/b1

op_basis = [op_one]
coeff_basis = [b1]

n_iter = 420
   

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

op_one = super_op(bh_hamil,op_zero)
op_one = op_one - trace_norm(op_zero,op_one)*op_zero
op_one = op_one - trace_norm(op_zero,op_one)*op_zero
b1 = norm(op_one)
op_one = op_one/b1
op_basis = [op_one]
coeff_basis = [b1]

for i in 2:n_iter
    new_op = super_op(bh_hamil,op_basis[i-1])
    
    new_op = new_op - trace_norm(op_zero,new_op)*op_zero
    for j in 1:i-1
        new_op = new_op - trace_norm(op_basis[j],new_op)*op_basis[j]
    end
    
    new_op = new_op - trace_norm(op_zero,new_op)*op_zero
    for j in 1:i-1
        new_op = new_op - trace_norm(op_basis[j],new_op)*op_basis[j]
    end
    
    new_coeff = norm(new_op)
    
    push!(op_basis,new_op/new_coeff)
    push!(coeff_basis,new_coeff)
end           



#y_axis = [coeff_basis[1:10]]
twofivecoup = coeff_basis
threecoup = coeff_basis
fourcoup = coeff_basis
twofiveop = op_basis
fourop = op_basis
inte = coeff_basis
x_axis = collect(1:n_iter)
plot(x_axis,[twofivecoup,fourcoup,inte],linewidth=1.3,label=[ "U/J=2.5" "U/J=4" "U/J=0"])
xlabel!("n")
ylabel!("Lanczos Coefficients")
savefig("/Users/priyoshankarpal/Desktop/coeff_BH_COMB.pdf")


Matrix(bh_hamil)
function time_evol_op(time, op, hamil)
    exp(1im*time*hamil)*op*exp(-1im*time*hamil)
end

times = range(0.0, stop=0.3, length=100)
time_axis = collect(times)
complexity = []
entropy =[]
for i in 1:length(time_axis)
    op_t = time_evol_op(time_axis[i], op_zero/norm(op_zero),Matrix(bh_hamil))
    sum = 0
    sum1 = 0
    phi_zero = trace_norm(op_zero,op_t)
    sum =sum + 0*abs(phi_zero)*abs(phi_zero)
    sum1 = sum1 - abs(phi_zero)*abs(phi_zero)*log(abs(phi_zero)*abs(phi_zero)) 
    for j in 1:n_iter
        phi_j = trace_norm(fourop[j],op_t)
        sum = sum + j*abs(phi_j)*abs(phi_j)
        sum1 = sum1 - abs(phi_j)*abs(phi_j)*log(abs(phi_j)*abs(phi_j))
    end
    push!(complexity,sum)
    push!(entropy, sum1)
end
c
complexity
plot(time_axis,[complexity, entropy],linewidth = 2.5,label =["complexity" "entropy"])
entropy
vline!([1.00], label="OTOC Decay", line=(color="black", linestyle=:dash))
xlabel!("t")
ylabel!("U/J=4")
savefig("/Users/priyoshankarpal/Desktop/complexity_BH_4early.pdf")



using LinearAlgebra
using Plots


# Compute eigenvalues
eigenvalues = eigvals(Matrix(bh_hamil))

# Sort the eigenvalues in ascending order
sorted_eigenvalues = sort(eigenvalues)

# Calculate level spacings
level_spacings = diff(sorted_eigenvalues)
histogram(level_spacings, bins=60)
num = length(level_spacings)
r_alpha = []
for i in 1:num
    if i!= num
        ra = level_spacings[i+1]/level_spacings[i]
        push!(r_alpha,min(ra,(1/ra)))
    end
end        
weak = r_alpha       
histogram([r_alpha], bins=20)
mean(r_alpha)

sort(level_spacings)
using Statistics
avg = mean(r_alpha)
# Plot the level spacing distribution
# Plot the level spacing distribution
histogram(r_alpha, bins=20, label="Level Spacing Distribution", legend=true)

# Customize the plot as needed
title!("Level Spacing Distribution for Quantum Chaos")
xlabel!("Level Spacing")
ylabel!("Frequency")

# Show the plot
plot()



    


