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



"""
Set Parameters

"""

N = 3  #number of sites in the chain  (int values only)                                      
M = 3   #number of bosons in the chain (int values only)
J = 4.0 #hopping paramter (float values only)
U = 8.5 # on-site potential (float values only) 
T =eltype(J)  #set data-type for the rest of the code
beta = 1.0  #inverse temperature

t_stop = 0.2
num_points = 30
np = pyimport("numpy")
times = np.linspace(0, t_stop, num_points)

#red_ham = [RBoseHubbard(N+1,M,J,U), RBoseHubbard(N,M,J,U)] 

"""
Step 1: Set the state of the system to contain only 1 boson
"""
pure_system = zeros(T, N+1, N+1)   #matrix of zeros for system dm that can contain from 0 to N bosons
pure_system[N, N] = 1.0            # The first element is 1, rest are zeros

"""
Step 2: Get the thermal density matrix 
"""
thermal_dm = thermal_state(beta,N, M, J, U)    #creates a block diagonal thermal dm with N-1 bosons and M-1 sites

"""
Step 3: Take the tensor product of the density matrix and thermal density matrix
"""
result_dm = sparse(kron(D, thermal_dm)) 

"""
Step 4: Limit the joint state to that of N bosons in total and apply the quench to set the intial state
(Quench process has been commented out for now. There are some errors when we perform quenching)
"""
#number_quench(N,M)
#number_quench_dag = SparseArrays.transpose(number_quench(N,M))
#result_dm
#limit_dm(result_dm,N,M)
#init_state = number_quench(N,M)*limit_dm(result_dm,N,M)*number_quench_dag

init_state = limit_dm(result_dm,N,M)
#init_state = init_state/tr(init_state)


"""
Define the required Hamiltonians
"""

H = BoseHubbard(N, M, J, U , :OBC)  #Full Hamiltonian of N bosons and M sites

sys_ham = RBoseHubbard(N, 2, 0.0, U)  #System Hamiltonian of N bosons and 1 site

bath_ham = RBoseHubbard.([N+1,N], M, J, U)  #Bath Hamiltonian of N bosons and M-1 sites
eigenvals, eigenvecs = eigen!(Matrix(bath_ham[2].H))  #eigenvalues and eigenvectors of the bath hamiltonian being used to trace over states





"""
Define empty lists for the bounds and otoc/renyi renyi_entropy
"""

renyi_ent_list =[]
bound_list=[]
bound_list_born=[]
spec_bound_list=[]
spec_bound_list_born=[]

@showprogress for (i,t) in enumerate(times)
    rho_t = time_evol_state(init_state, H, t )
    rho_int_t = interaction_picture(NBasis(N,M), U, rho_t)
    
    rho_S = partial_trace_bath(rho_int_t, N, M)
    push!(renyi_ent_list, renyi_entropy(rho_S))

    rho_B = partial_trace_system(rho_int_t, size(init_state,1),N,M)

   

    qsl = QSL_OTOC(t, J, N, sys_ham, bath_ham, eigenvecs, rho_B)
    qsl_born = QSL_OTOC(t, J, N, sys_ham, bath_ham, eigenvecs, thermal_dm)

    push!(bound_list, real(qsl.state_bound))
    push!(spec_bound_list, real(qsl.spectral_bound))
    push!(bound_list_born, real(qsl_born.state_bound))
    push!(spec_bound_list_born, real(qsl_born.spectral_bound))
    

end   

plot(times, [exp.(-2*real(bound_list)),exp.(-2*real(bound_list_born)), exp.(-spec_bound_list) , exp.(-spec_bound_list_born),exp.(-renyi_ent_list)] )


