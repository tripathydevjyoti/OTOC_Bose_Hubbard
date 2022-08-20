@testset "Basis" begin
    M = N = 3

    B = Basis(N, M)
    NB = NBasis(N, M)

    SB = [
        [3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0],
        [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1],
        [0, 1, 2], [0, 0, 3]
    ]
    I = B.sub_basis_indices

    @test length(B.eig_vecs) == length(B.tags) == (M + 1) ^ N == B.dim
    @test eltype(B.eig_vecs) == Vector{Int} == eltype(NB.eig_vecs)
    @test all(length.(B.eig_vecs) .== M) && all(length.(NB.eig_vecs) .== M)
    @test length(I) == sub_basis_dim(N, M) == length(NB.eig_vecs) == NB.dim
    @test Set(B.eig_vecs[I]) == Set(SB) == Set(NB.eig_vecs)
    @test all(sum.(B.eig_vecs[I]) .== N) && all(sum.(NB.eig_vecs) .== N)
    @test B.N == N == NB.N && B.M == M == NB.M
    @test issorted(B.tags) && issorted(NB.tags)
end
