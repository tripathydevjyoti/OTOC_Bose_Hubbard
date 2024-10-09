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
using UnicodePlots
using ExponentialUtilities




N=4
M=4
U=9.0
J = 4.0 #hopping paramter (float values only)
beta = 1.0

T =eltype(J)
struct RBoseHubbard{T <: Number}
    basis::Union{Basis, NBasis, Vector{NBasis}, RBasis}
    lattice::LabelledGraph
    H::SparseMatrixCSC{T, Int64}
end

# Outer constructor to infer T automatically
function RBoseHubbard(N::Int, M::Int, J::T, U::T) where T <: Number
    RH =BlockDiagonal([BoseHubbard(i, M-1, J, U, :OBC).H for i in N:-1:0] ) # Vector of Hamiltonians
    RB = RBasis(N, M)  # Basis initialization
    lattice = chain(M, J, U, :OBC)  # Lattice initialization
    return RBoseHubbard{T}(RB, lattice, RH)  # Return an instance of RBoseHubbard
end

red_ham = RBoseHubbard(N,M,J,U)







#init_state = eigenvecs[:,1]

function thermal_dm(ham::RBoseHubbard{T}, beta::T) where T <: Number
   

    eigenvals, eigenvecs = eigen(Matrix(ham.H))
    
    exp_eigenvals = exp.(-beta * eigenvals)

    Z = sum(exp_eigenvals)
    
    ρ = Diagonal(exp_eigenvals / Z)
    return round.(eigenvecs*ρ*eigenvecs' , digits=10)
    

end













function partial_trace(init_dm, subsys_size)
    
    r_dm = zeros(ComplexF64, subsys_size, subsys_size )
    for (i,_) in enumerate(1:subsys_size)
        for (j,_) in enumerate(1:subsys_size)
            if init_dm[i,j] != 0
                if NBasis(N,M).eig_vecs[i][1] == NBasis(N,M).eig_vecs[j][1]
      
                    r_vec_i = deleteat!(NBasis(N,M).eig_vecs[i], 1)
                    r_vec_j = deleteat!(NBasis(N,M).eig_vecs[j], 1)
    
                    index1 = findfirst(x -> x ==r_vec_i, RBasis(N,M).eig_vecs)
                    index2 = findfirst(x -> x ==r_vec_j, RBasis(N,M).eig_vecs)
              
                    r_dm[index1, index2] = r_dm[index1, index2] + init_dm[i,j]
                end    
 
             end
         end        
    end 
    return r_dm 
end    


np = pyimport("numpy")
t_stop = 40.0
num_points = 200
times = np.linspace(0, t_stop, num_points)

function time_evolution(t, H, state)

    time = -1im*t
    time_evol_state = expv(time, H, State(state, H.basis) )
end    




function two_time_corr(
     H::RBoseHubbard{S},eigss,  time::T, rho;  kwargs=()
    ) where{S, T <:Real}
    
    trsum1 = 0
    trsum2 = 0
    
   
    for i in 1:length(H.basis.dim)
        
       
        
        gamma1 = bath_exact(time, time, H, 1, 1, State(eigss[:,i] ,H.basis), rho )
        trsum1 = trsum1 + gamma1 

        gamma2 = bath2_exact(time, time, H, 1, 1, State(eigss[:,i], H.basis), rho)
        
        trsum2 = trsum2 + gamma2
       
    end
    
    corr_arr = [trsum1, trsum2 ] 
return corr_arr
end

get_index(RBasis(3,3),[2,1])


@showprogress for U in [9.0]
        
    #BH hamiltonian 
    eigenvals, eigenvecs = eigen(Matrix(red_ham.H))
    init_dm =thermal_dm(red_ham, beta)
    
    
    subsys_size = red_ham.basis.dim
    
    twopt1 =[] #array to store values for Γ1
    twopt2 =[] #array to store values for Γ2
  

    #eigenvals1, eigenvecs1 = eigen(Matrix(H[2].H))
    for (_, t) in enumerate(times)
        τ = 1im*t
        time_evol_dm = exponential!(τ*Matrix(red_ham.H))*init_dm*exponential!(-τ*Matrix(red_ham.H))
        red_dm = partial_trace(time_evol_dm, subsys_size)
        arr = two_time_corr( red_ham, eigenvecs, t, red_dm)
        push!(twopt1, arr[1])
        push!(twopt2, arr[2])
    end 
    print(U)
    
    
    filename1 = @sprintf("N_%d_L_%d_BH_0temp_U_%.1f_J_%.1f_t_%.1f_num_points_%.1f_Gamma1.npy", N, M, U, J, t_stop, num_points)
    filename2 = @sprintf("N_%d_L_%d_BH_0temp_U_%.1f_J_%.1f_t_%.1f_num_points_%.1f_Gamma2.npy", N, M, U, J, t_stop, num_points)
    np.save(filename1,twopt1)
    np.save(filename2,twopt2)
    
end    

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