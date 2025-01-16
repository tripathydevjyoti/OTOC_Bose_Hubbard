export
   thermal_state,
   partial_trace_system,
   partial_trace_bath,
   two_time_corr,
   thermal_corr

"""
$(TYPEDSIGNATURES)
"""   






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




function partial_trace_system(init_dm, subsys_size,N,M)
    
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

function thermal_corr(
    beta::T, H::Vector{BoseHubbard{S}}, time::T, partition::T, eigenvals, eigenvecs;  kwargs=()
    ) where{S, T <:Real}
    
    trsum1 = 0
    trsum2 = 0
    
   
    for (i, energy) in enumerate(eigenvals)
        
        gibbs = exp(-beta*energy)
        
        gamma1 = gibbs*bath(time, 0.0, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        trsum1 = trsum1 + gamma1 

        gamma2 = gibbs*bath2(time,0.0, H, 1, 1, State(eigenvecs[:, i], H[2].basis) )
        trsum2 = trsum2 + gamma2
       
    end
    
    corr_arr = [trsum1/ partition, trsum2/partition ] 
return corr_arr
end


function two_time_corr(
    H::Vector{RBoseHubbard{S}},eigss,  time::Vector{T}, rho;  kwargs=()
   ) where{S, T <:Real}
   
   trsum1 = 0
   trsum2 = 0
   
  
   for i in 1:length(H[2].basis.dim)
       
      
       
       gamma1 = bath_exact(time[1], time[2], H, 1, 1, State(eigss[:,i] ,H[2].basis), rho )
       trsum1 = trsum1 + gamma1 

       gamma2 = bath2_exact(time[1], time[2], H, 1, 1, State(eigss[:,i], H[2].basis), rho)
       
       trsum2 = trsum2 + gamma2
      
   end
   
   corr_arr = [trsum1, trsum2 ] 
return corr_arr
end

"""
function two_time_corr(
H::Vector{RBoseHubbard{S}},eigss,  time::T, rho;  kwargs=()
  ) where{S, T <:Real}
  two_time = [time, 0]
  two_time_corr(H, eigss, two_time,rho)
end   
"""