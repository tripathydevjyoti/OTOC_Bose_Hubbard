M = N = 3

@testset "Basis (conserved particles)" begin
    dim = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))
    B_ref = [
        [3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0],
        [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1],
        [0, 1, 2], [0, 0, 3]
    ]
    B = Basis(M, N; constraint=:conserved_particles)

    @test length(B.eig_vecs) == length(B.tags) == dim
    @test eltype(B.eig_vecs) == Vector{Int}
    @test all(length.(B.eig_vecs) .== M)
    @test Set(B.eig_vecs) == Set(B_ref)
    @test all(sum.(B.eig_vecs).== N)
    @test B.N == N
    @test B.M == M
end

@testset "Basis (full)" begin
    B = Basis(M, N)
    dim = (N + 1) ^ M

    @test length(B.eig_vecs) == length(B.tags) == dim
    @test eltype(B.eig_vecs) == Vector{Int}
    @test B.N == N
    @test B.M == M
end
