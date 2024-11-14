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
U=8.0
T =eltype(J)
np = pyimport("numpy")

t_stop = 40.0
num_points = 200
times = np.linspace(0, t_stop, num_points)
red_ham = [RBoseHubbard(N+1,M,J,U), RBoseHubbard(N,M,J,U)] 

tensor_basis(N,M)

function thermal_state(beta::T, ham::RBoseHubbard{T}) where T<:Real
    thermal_mat = sparse(exponential!(Matrix(-beta*ham.H)))
    part_func = tr(thermal_mat)

    return thermal_mat/part_func
end


sys_basis = RBasis(N,2).eig_vecs
bath_basis = RBasis(N,M).eig_vecs
products  = collect.(Iterators.product(sys_basis,bath_basis))

states = vec(np.array([vcat(p...) for p in products].tranpose())

tensor_basis(N,M)
       






 



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