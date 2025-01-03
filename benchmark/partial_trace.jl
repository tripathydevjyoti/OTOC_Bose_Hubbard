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
using ProgressMeter
using Printf
using SparseArrays
using BlockDiagonals
using ExponentialUtilities
using LogExpFunctions





N=3
M=3
J = 4.0 #hopping paramter (float values only)
beta = 1.0
U = 9.0
T =eltype(J)
np = pyimport("numpy")

t_stop = 0.25
num_points = 30
times = np.linspace(0, t_stop, num_points)
red_ham = [RBoseHubbard(N+1,M,J,U), RBoseHubbard(N,M,J,U)] 

hamil = BoseHubbard(N,M,J,U,:OBC).H
eigvals, eigvecs = eigen!(Matrix(hamil))
diff(sort(eigvals))
function thermal_state(beta::T, N::Int, M::Int, J::T, U::T) where T<:Real
    # Original Hamiltonian
    ham = BoseHubbard(N-1, M-1, J, U, :OBC).H

# Compute the thermal density matrix for the Hamiltonian
    thermal_mat = sparse(exponential!(Matrix(-beta * ham)))
    part_func = tr(thermal_mat)
    thermal_dm = thermal_mat / part_func  # Normalize to get the density matrix

# First block: Zero matrix of size NBasis(N, M-1)
    first_block_dim = NBasis(N, M-1).dim
    first_block = spzeros(first_block_dim, first_block_dim)  # Sparse zero matrix

# Second block: The Hamiltonian
    second_block = thermal_dm

# Remaining zero blocks
    zero_blocks = []
    for i in (N-2):-1:0
        dim = NBasis(i, M-1).dim  # Dimension of the current zero block
        push!(zero_blocks, spzeros(dim, dim))  # Add sparse zero matrix of appropriate size
    end

    
# Combine all blocks into a block diagonal matrix
    blocks = [first_block, second_block, zero_blocks...]  # Concatenate all blocks
    block_diag_matrix = BlockDiagonal(blocks)

# Return the block diagonal matrix
return block_diag_matrix

end


#print(thermal_state(beta, N, M, J, U))




function limit_dm(rho::SparseMatrixCSC{T, Int64}, N::Int, M::Int) where T<:Number
    spinfo = findnz(rho)

    dims = NBasis(N,M).dim
    limit_rho = zeros(dims, dims)
    for (i,val) in enumerate(findnz(result_dm)[3])
        
        index1 = get_index(NBasis(N,M), tensor_basis(N,M).eig_vecs[spinfo[1][i]])
        index2 = get_index(NBasis(N,M), tensor_basis(N,M).eig_vecs[spinfo[2][i]])
        
        limit_rho[index1, index2] = val
    end
    
    
    return sparse(limit_rho)
end    



function number_quench(N,M)
    size = NBasis(N,M).dim
    vecs = NBasis(N,M).eig_vecs
    quench = zeros(size,size)
    for (i, _) in enumerate(vecs)
        for (j, state) in enumerate(vecs)
            quench[i,j] = state[1]
        end    
    end

    
    sparse(quench)
end

function time_evol_state(rho::SparseMatrixCSC{T, Int64}, bh::BoseHubbard{T}, time::T) where T<: Number
    τ = -1im*time
    U = exponential!(Matrix(τ*bh.H))
    U_dag = adjoint(U)

    U*rho*U_dag

end    

function partial_trace_bath(init_dm, N, M)
    
    r_dm = zeros(ComplexF64, N+1, N+1 )
    for (i,_) in enumerate(1:size(init_dm, 1))
        for (j,_) in enumerate(1:size(init_dm, 2))
            if init_dm[i,j] != 0
                if NBasis(N,M).eig_vecs[i][2:end] == NBasis(N,M).eig_vecs[j][2:end]
    
                    index1 = NBasis(N,M).eig_vecs[i][1]+1
                    index2 = NBasis(N,M).eig_vecs[j][1]+1
                    r_dm[index1, index2] = r_dm[index1, index2] + init_dm[i,j]
                end    
 
             end
         end        
    end 
    return sparse(r_dm) 
end 

function renyi_entropy(rho)
    rho_sq = rho*rho
    tr_rho_sq = real(tr(rho_sq))
    
    return -log(tr_rho_sq)

end   



D = zeros(T, N+1, N+1)
D[N, N] = 1.0  # The first element is 1, rest are zeros

    # Step 2: Get the thermal density matrix (assuming thermal_state is defined)
thermal_dm = thermal_state(beta,N, M, J, U)

    # Step 3: Take the tensor product of the density matrix and thermal density matrix
result_dm = sparse(kron(D, thermal_dm)) 

number_quench(N,M)
number_quench_dag = SparseArrays.transpose(number_quench(N,M))
result_dm
limit_dm(result_dm,N,M)
#init_state = number_quench(N,M)*limit_dm(result_dm,N,M)*number_quench_dag
init_state = limit_dm(result_dm,N,M)
#init_state = init_state/tr(init_state)
print(tr(init_state))
H = BoseHubbard(N, M, J, U , :OBC)

rho_t = time_evol_state(init_state, H, 2.0 )

renyi_ent_list =[]
renyi_ent_list2=[]
for (_,t) in enumerate(times)
    rho_t = time_evol_state(init_state, H, t )
    rho_B = partial_trace(rho_t, size(init_state,1),N,M)
    rho_S = partial_trace_bath(rho_t, N, M)
    println(size(rho_B))
    push!(renyi_ent_list, renyi_entropy(rho_S))
    push!(renyi_ent_list2, renyi_entropy(rho_B))
end    


bath_ham = RBoseHubbard.([N+1,N], M, J, U)
eigenvals, eigenvecs = eigen!(Matrix(bath_ham[2].H))

two_time_corr(bath_ham, eigenvecs , [2.0,0.0], thermal_dm)

function create_annihilation_creation_descending(N::Int)
    # Initialize (N+1)x(N+1) matrices
    a = zeros(ComplexF64, N+1, N+1)  # Annihilation operator
    adag = zeros(ComplexF64, N+1, N+1)  # Creation operator

    # Populate the matrices
    for n in 1:N
        a[n+1, n] = sqrt(N - n + 1)  # a lowers |N-n+1⟩ to |N-n⟩
        adag[n, n+1] = sqrt(N - n + 1)  # a† raises |N-n⟩ to |N-n+1⟩
    end

    return a, adag
end


a, adag = create_annihilation_creation_descending(N)

function time_evol_jump(time,op)
    τ =-1im*time
    prop = exponential!(Matrix(τ*RBoseHubbard(N, 2, 0.0, U).H))
    prop_dag = adjoint(prop)
    evol_op = prop_dag*op*prop
    norm(evol_op, Inf)
end  

function integrand(time1, time2, J)
    bath =two_time_corr(bath_ham, eigenvecs , [time1, time2], thermal_dm)
    norm1 = time_evol_jump(time1, adag)*time_evol_jump(time2, a)
    norm2 = time_evol_jump(time1, a)*time_evol_jump(time2, adag)
    term1 = (bath[1]+conj(bath[2]))*2*norm2
    term2 = (conj(bath[1])+bath[2])*2*norm1

    J*J*(term1+term2)
end 

function double_integral(t, J)
    inner_integral(time1) = quadgk(time2 -> integrand(time1, time2, J), 0, time1)[1]
    quadgk(inner_integral, 0, t)[1]
end

bound_list=[]
for (_,t) in enumerate(times)
    push!(bound_list, real(double_integral(t,J)))
    print(t)

end   
bound_list
plot(times, [exp.(-real(bound_list)), exp.(-renyi_ent_list)])

plot(times, [real(bound_list),renyi_ent_list ])

println("Annihilation Operator (a) in descending basis:")
println(a)
println("\nCreation Operator (a†) in descending basis:")
println(adag)

renyi_entropy(D)
renyi_ent_list2
renyi_ent_list
plot(times, [renyi_ent_list, renyi_ent_list2])
plot(times, exp.(-renyi_ent_list))

eigen(Matrix(init_state)).values





"""
sys_basis = RBasis(N,2).eig_vecs
bath_basis = RBasis(N,M).eig_vecs
products  = collect.(Iterators.product(sys_basis,bath_basis))

states = vec(np.array([vcat(p...) for p in products].tranpose())

tensor_basis(N,M)
"""     






 



"""
using ExponentialUtilities


function time_evol_rho(t, rho, H)
    τ = -1im*t
    mat = Matrix(h.H)
    U_t = exponential!(τ*mat)
    U_t_dag = conj(U_t)
    return U_t*rho*U_t_dag

end

red_ham = [RBoseHubbard(N+1,M,J,U), RBoseHubbard(N,M,J,U)] 
time_evol_rho(1.0, thermal_dm, h)
@showprogress for U in [4.0]
        
    h = BoseHubbard(N,M,J,U,:OBC)
    red_ham = [RBoseHubbard(N+1,M,J,U), RBoseHubbard(N,M,J,U)] 
    rho_beta = exponential!(-beta*Matrix(h.H))
    z = tr(exponential!(-beta*Matrix(h.H)))
    thermal_dm = rho_beta/z
    
    subsys_size = h.basis.dim
    
    twopt1 =[] #array to store values for Γ1
    twopt2 =[] #array to store values for Γ2
  
    for (_, t) in enumerate(times)
        
        time_evol_dm = time_evol_rho(t,thermal_dm,h)
        red_dm = partial_trace(time_evol_dm, subsys_size,N,M)
        arr = two_time_corr( red_ham, eigenvecs, t, red_dm)
        push!(twopt1, arr[1])
        push!(twopt2, arr[2])
    end 
    print(U)
    
    
    #filename1 = @sprintf("N_%d_L_%d_BH_0temp_U_%.1f_J_%.1f_t_%.1f_num_points_%.1f_Gamma1.npy", N, M, U, J, t_stop, num_points)
    #filename2 = @sprintf("N_%d_L_%d_BH_0temp_U_%.1f_J_%.1f_t_%.1f_num_points_%.1f_Gamma2.npy", N, M, U, J, t_stop, num_points)
    #np.save(filename1,twopt1)
    #np.save(filename2,twopt2)
    
end    

"""

#twopt1 = np.load("N_5_L_5_BH_0temp_U_8.0_J_4.0_t_40.0_num_points_200.0_Gamma1.npy")
#using Plots
#plot(times, real(twopt1))




"""
using Plots
plot(times, real(twopt1))

H = BoseHubbard.([N+1, N, N-1,N-2], M, J, 4.0, :OBC)
eigvals1, eigvecs1 = eigen(Matrix(H[2].H))
for (_,vec) in enumerate(eigvecs1)
    println(vec)
end    
length(eigvecs1)
State(eigvecs1[:,1],R_BH2.basis)
H[2].basis.eig_vecs[1]
State(H[2].basis.eig_vecs[1], H[2].basis)
State([1.0], [H[2].basis.eig_vecs[1]])
H = BoseHubbard.([N+1, N, N-1,N-2], M, J, 2.0, :OBC) #BH hamiltonian 
for (i, vec) in enumerate(H[2].basis.eig_vecs)
    print(eltt)
end    

H = BoseHubbard.([N+1, N, N-1,N-2], M, J, U, :OBC) #BH hamiltonian 
eigenvals, eigenvecs = eigen(Matrix(H[2].H))
init_state = eigenvecs[:,1]
    #init_state =dense( State([1.0],[[1,1,1,1,1,1]]), H[2].basis)
subsys_size = H[2].basis.dim
init_dm = partial_trace(init_state * init_state', subsys_size)
twopt1 =[] #array to store values for Γ1
twopt2 =[] #array to store values for Γ2
R_BH1 = BoseHubbard(RBasis(N+1,M-1), chain(M-1, J, U, :OBC))
R_BH2 = BoseHubbard(RBasis(N,M-1), chain(M-1, J, U, :OBC))
R_BH3 = BoseHubbard(RBasis(N-1,M-1), chain(M-1, J, U, :OBC))
time1 = 0.1
time2 = 0.1
τ = -1im * time1
s = -1im * time2
H1 = [R_BH1, R_BH2, R_BH3]

    #therm_ket = expv(-beta, H[2], state
state = State([1.0],[H1[2].basis.eig_vecs[1]])

evol_bra = expv((-τ), H1[2], state)
bdag_bra = expv( τ , H1[1], create(State(evol_bra,H1[2].basis), 1))
state1= State(bdag_bra, H1[1].basis)
n = length(state1.eig_vecs)
vecs = Vector(undef, n)
Threads.@threads for k ∈ 1:n
    ket = state1.eig_vecs[k]
    println(ket)
    vecs[k] = ket[1] > 0 ? operate(ket, 1, -1) : 0
    println(vecs[k])
end
vecs
K = findall(!iszero, vecs)
kets = [0.0+0.0im, 0.0+0.0im]
State(state1.coeff[K] .* sqrt.(getindex.(kets, 1) .+ 1), kets)


half_bra = dense(destroy(State(bdag_bra, H1[1].basis), 1), H1[2].basis)
""" 




"""
#init_dm = fill(1/3, (length(eigenvals), length(eigenvals)))
r_dm = zeros(ComplexF64, length(eigenvals), length(eigenvals) )
for (i,_) in enumerate(1:length(eigenvals))
    for (j,_) in enumerate(1:length(eigenvals))
        if init_dm[i,j] != 0
            if NBasis(N,M).eig_vecs[i][1] == NBasis(N,M).eig_vecs[j][1]
      
                r_vec_i = deleteat!(NBasis(N,M).eig_vecs[i], 1)
                r_vec_j = deleteat!(NBasis(N,M).eig_vecs[j], 1)
    
                index1 = findfirst(x -> x ==r_vec_i, RBasis(N,M-1).eig_vecs)
                index2 = findfirst(x -> x ==r_vec_j, RBasis(N,M-1).eig_vecs)
              
                r_dm[index1, index2] = r_dm[index1, index2] + init_dm[i,j]
            end    
 
        end
    end        
end    

print(Matrix((r_dm)))
"""