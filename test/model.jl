using LinearAlgebra
using SparseArrays

commutator(A, B) = A * B .- B * A

M = N = 3
D = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))

B_ref = [
    [3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0],
    [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1],
    [0, 1, 2], [0, 0, 3]
]
B = Basis(M, N)

@testset "Basis" begin
    @test length(B.eig_vecs) == length(B.tags) == D
    @test eltype(B.eig_vecs) == Vector{Int}
    @test all(length.(B.eig_vecs) .== M)
    @test Set(B.eig_vecs) == Set(B_ref)
    @test all(sum.(B.eig_vecs).== N)
    @test B.N == N
    @test B.M == M
end

@testset "Operate (destroy and create)" begin
    @test operate([3, 0, 0], 1, :destroy) == [2, 0 ,0]
    @test operate([3, 0, 0], 2, :create) == [3, 1 ,0]
    @test destroy_and_create([2, 0, 1], 1, 2) == [1, 1, 1]
end

@testset "Dense occupation operators" begin
    for T ∈ (Float16, Float32, Float64)
        n = occupation(T, B, :dense)
        @test eltype(n) == T
        @test isdiag(n)
        @test transpose(n) == n
        @test all(diag(n) .== B.N)
        @test n == sum(occupation(B, i, :dense) for i ∈ 1:M)
    end
end

@testset "Sparse occupation operators" begin
    for T ∈ (Float16, Float32, Float64)
        for i ∈ 1:M
            n_i = occupation(T, B, i, :sparse)
            @test eltype(n_i) == T
            @test isdiag(n_i)
            @test transpose(n_i) == n_i
        end
    end
end

@testset "Dense Hamiltonian" begin
    T = Float64
    J = T(1)
    U = T(1/2)
    #H = hamiltonian(M, N, J, U, :OBC)
    #@test size(H) = (D, D)
end
