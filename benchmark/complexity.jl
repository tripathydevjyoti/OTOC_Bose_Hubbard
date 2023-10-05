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


for i in 1:dim_ham
    op_zero[i,i] = float(H[1].basis.eig_vecs[i][site])
end    


function super_op(H,op)
    return H*op - op*H

end  

op_one = super_op(H[1].H,op_zero)

function norm(op)
    sqrt(tr(adjoint(op)*op))
end 





b1 = norm(op_one)
op_one = op_one/b1

op_basis = [op_one]
coeff_basis = [b1]

n_iter = 10

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

x_axis = [1,2,3,4,5,6,7,8,9,10]

plot(x_axis, coeff_basis)

    