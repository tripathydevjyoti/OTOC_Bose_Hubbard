export
   thermal_dm,
   partial_trace,
   two_time_corr,
   thermal_corr

"""
$(TYPEDSIGNATURES)
"""   






function thermal_dm(ham::RBoseHubbard{T}, beta::T) where T <: Number
   

    eigenvals, eigenvecs = eigen(Matrix(ham.H))
    
    exp_eigenvals = exp.(-beta * eigenvals)

    Z = sum(exp_eigenvals)
    
    ρ = Diagonal(exp_eigenvals / Z)
    return round.(eigenvecs*ρ*eigenvecs' , digits=10)
    

end




function partial_trace(init_dm, subsys_size,N,M)
    
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