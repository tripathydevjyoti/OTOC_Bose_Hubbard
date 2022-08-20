M = N = 3
B = Basis(M, N)
I = B.sub_basis_indices
D = sub_basis_dim(N, M)

@testset "Sparse occupation operators" begin
    for T ∈ (Float16, Float32, Float64)
        for i ∈ 1:M
            n_i = occupation(T, B.eig_vecs, i)
            sub_n_i = occupation(T, B.eig_vecs[I], i)
            @test length(nonzeros(n_i)) == B.dim == nnz(n_i)
            @test length(nonzeros(sub_n_i)) == D == nnz(sub_n_i)
            @test sub_n_i == n_i[I, I]
            @test eltype(n_i) == eltype(sub_n_i) == T
            @test isdiag(n_i) && isdiag(sub_n_i)
            @test transpose(n_i) == n_i && transpose(sub_n_i) == sub_n_i
            @test size(n_i) == (B.dim, B.dim)
            @test size(sub_n_i) == (D, D)
        end
        @test occupation(T, B.eig_vecs) == sum(occupation(B.eig_vecs, i) for i ∈ 1:M)
        @test occupation(T, B.eig_vecs[I]) == sum(occupation(B.eig_vecs[I], i) for i ∈ 1:M)
    end
end
