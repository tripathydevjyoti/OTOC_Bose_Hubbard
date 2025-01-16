export
   QSL_OTOC,
   time_evol_state,
   interaction_picture,
   renyi_entropy,
   time_evol_jump,
   time_evol_jump_norm,
   limit_dm,
   number_quench,
   create_annihilation_creation_descending

"""
$(TYPEDSIGNATURES)
""" 
   
  


struct QSL_OTOC{T <: Real}
    state_bound::T
    spectral_bound::T

    function QSL_OTOC(t::T, J::T, N::Int, sys_ham::RBoseHubbard{T}, bath_ham::Vector{RBoseHubbard{T}}, eigenvecs::Matrix{T},rhob::Union{Matrix{ComplexF64}, BlockDiagonal{Float64, SparseMatrixCSC{Float64, Int64}}}) where T <: Real
        state_bound = real(compute_state_bound(t, J, N, sys_ham, bath_ham, eigenvecs, rhob))
        spectral_bound = real(compute_spectral_bound(t, J, N, sys_ham, bath_ham, eigenvecs, rhob))
        new{T}(state_bound, spectral_bound)
    end

    
    function compute_state_bound(t::T, J::T, N::Int, sys_ham::RBoseHubbard{T}, bath_ham::Vector{RBoseHubbard{T}}, eigenvecs::Matrix{T}, rhob::Union{Matrix{ComplexF64}, BlockDiagonal{Float64, SparseMatrixCSC{Float64, Int64}}}) where T <: Real
        function integrand(time1::T, time2::T, J::T) where T <: Real
            bath = two_time_corr(bath_ham, eigenvecs, [time1, time2], rhob)
            a, adag = create_annihilation_creation_descending(N)
            norm1 = time_evol_jump_norm(time1, adag, sys_ham) * time_evol_jump_norm(time2, a, sys_ham)
            norm2 = time_evol_jump_norm(time1, a, sys_ham) * time_evol_jump_norm(time2, adag, sys_ham)
            term1 = (bath[1] + conj(bath[2])) * 2 * norm2
            term2 = (conj(bath[1]) + bath[2]) * 2 * norm1
            J * J * (term1 + term2)
        end

        inner_integral(time1::T) = quadgk(time2 -> integrand(time1, time2, J), 0, time1)[1]
        quadgk(inner_integral, 0, t)[1]
    end

   
    function compute_spectral_bound(t::T, J::T, N::Int, sys_ham::RBoseHubbard{T}, bath_ham::Vector{RBoseHubbard{T}}, eigenvecs::Matrix{T}, rhob::Union{Matrix{ComplexF64}, BlockDiagonal{Float64, SparseMatrixCSC{Float64, Int64}}}) where T <: Real
        function liouv_superop(time1::T, time2::T, rhob::Union{Matrix{ComplexF64}, BlockDiagonal{Float64, SparseMatrixCSC{Float64, Int64}}}) where T <: Real
            a, adag = create_annihilation_creation_descending(N)
            a_t1 = time_evol_jump(time1, a, sys_ham)
            a_t2 = time_evol_jump(time2, a, sys_ham)
            adag_t1 = time_evol_jump(time1, adag, sys_ham)
            adag_t2 = time_evol_jump(time2, adag, sys_ham)

            bath_corr = two_time_corr(bath_ham, eigenvecs, [time1, time2], rhob)
            iden = zeros(N + 1, N + 1)
            superop1 = (J * J) * (bath_corr[1]) * (kron(adag_t2, a_t1) - kron(iden, adag_t2 * a_t1))
            superop2 = (J * J) * (conj(bath_corr[1])) * (kron(adag_t1, a_t2) - kron(adag_t1 * a_t2, iden))
            superop3 = (J * J) * (bath_corr[2]) * (kron(a_t2, adag_t1) - kron(iden, a_t2 * adag_t1))
            superop4 = (J * J) * (conj(bath_corr[2])) * (kron(a_t1, adag_t2) - kron(a_t1 * adag_t2, iden))

            superop = superop1 + superop2 + superop3 + superop4
            superop - adjoint(superop)
        end

        function inner_integral_spectral(time1::T) where T <: Real
            norm(quadgk(time2 -> liouv_superop(time1, time2, rhob), 0, time1)[1], Inf)
        end

        quadgk(time1 -> inner_integral_spectral(time1), 0, t)[1]
    end
end



function time_evol_state(rho::SparseMatrixCSC{T, Int64}, bh::BoseHubbard{T}, time::T) where T<: Number
    τ = -1im*time
    U = exponential!(Matrix(τ*bh.H))
    U_dag = adjoint(U)
   
    
    U*rho*U_dag
end  

function interaction_picture(basis, U, rho)
    dims = basis.dim
    rho_int_pic = zeros(ComplexF64, dims, dims)
    for (i,_) in enumerate(1:dims)
        for (j,_) in enumerate(1:dims)
            n1 = basis.eig_vecs[i][1]
            n2 = basis.eig_vecs[j][1]
            rho_int_pic[i,j] = rho[i,j]*exp(1im*0.5*U*n1*(n1-1))*exp(-1im*0.5*U*n2*(n2-1))
        end
    end
    rho_int_pic
end


function renyi_entropy(rho)
    rho_sq = rho*rho
    tr_rho_sq = real(tr(rho_sq))
    
    return -log(tr_rho_sq)

end   

function time_evol_jump_norm(time::T,op, ham::RBoseHubbard{T}) where T<:Real
    τ =-1im*time
    prop = exponential!(Matrix(τ*ham.H))
    prop_dag = adjoint(prop)
    evol_op = prop_dag*op*prop
    norm(evol_op, Inf)
end  

function time_evol_jump(time::T,op, ham::RBoseHubbard{T}) where T<:Real
    τ =-1im*time
    prop = exponential!(Matrix(τ*ham.H))
    prop_dag = adjoint(prop)
    evol_op = prop_dag*op*prop
   
end  


function limit_dm(rho::SparseMatrixCSC{T, Int64}, N::Int, M::Int) where T<:Number
    spinfo = findnz(rho)

    dims = NBasis(N,M).dim
    limit_rho = zeros(dims, dims)
    for (i,val) in enumerate(findnz(rho)[3])
        
        index1 = get_index(NBasis(N,M), tensor_basis(N,M).eig_vecs[spinfo[1][i]])
        index2 = get_index(NBasis(N,M), tensor_basis(N,M).eig_vecs[spinfo[2][i]])
        
        limit_rho[index1, index2] = val
    end
    
    
    return sparse(limit_rho)
end  

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

function number_quench(N,M)
    size = NBasis(N,M).dim
    vecs = NBasis(N,M).eig_vecs
    quench = zeros(size,size)
    for (i, state) in enumerate(vecs)
        
        for (j,_) in enumerate(vecs)
            quench[i,j] = state[1]
        end    
    end
    trace = tr(quench)

    
    sparse(quench/trace)
end
