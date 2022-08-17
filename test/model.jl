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
        @test occupation(T, B) == occupation(T, B, :sparse) |> Array
        for i ∈ 1:M
            n_i = occupation(T, B, i, :sparse)
            @test eltype(n_i) == T
            @test isdiag(n_i)
            @test transpose(n_i) == n_i
            @test n_i == occupation(T, B, i, :sparse) |> Array
        end
    end
end

#=
@testset "Dense creation and annihilation operators" begin
    for i ∈ 1:M, j ∈ 1:M
        a_i = annihilation(B, i, :dense)
        ap_j = creation(B, j, :dense)

        @test iszero(tr(a_i)) && iszero(tr(a_i))
        @test size(a_i) == size(ap_j) == (D, D)

        # [a_i^†, a_j^†] == [a_i, a_j] == 0
        @test iszero(commutator(a_i, ap_j |> transpose))
        @test iszero(commutator(a_i |> transpose, ap_j))

        # [a_i, a_j^†] == δ_ij
        if i != j
            @test iszero(commutator(a_i, ap_j))
        else
            #@test commutator(a_i, ap_j) == Matrix(1.0 * I, D, D)
            println(iszero(ap_j))
            println(iszero(a_i))
            n_i = occupation(B, i, :dense)
            @test transpose(a_i) == ap_j
            @test n_i == ap_j * a_i
        end
    end
end

@testset "Dense creation and annihilation operators" begin
    for i ∈ 1:M, j ∈ 1:M
        a_i = annihilation(B, i)
        ap_j = creation(B, j)

        @test issparse(a_i)
        @test issparse(ap_j)

        @test a_i == annihilation(B, i, :sparse)
        @test ap_j == creation(B, j, :sparse)

        @test tr(a_i) == tr(ap_j) == 0
        @test size(a_i) == size(ap_j) == (D, D)
        if i == j @test transpose(a_i) == ap_j end
    end
end
=#
