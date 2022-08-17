using LinearAlgebra
using SparseArrays

δ(i::Int, j::Int, D::Int) = i == j ? Matrix(1.0 * I, D, D) : zeros(D, D)
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

@testset "Sparse creation and annihilation operators" begin
    for i ∈ 1:2, j ∈ 1:1
        a_i = creation(B, i)
        ap_j = annihilation(B, j)

        @test issparse(a_i)
        @test issparse(ap_j)
        @test size(a_i) == size(ap_j) == (D, D)
        if i == j @test transpose(a_i) == ap_j end

        # [a_i, a_j] == 0
        b_i = a_i |> Array
        b_j = ap_j |> Array |> transpose

        @test commutator(b_i, b_j) ≈ zeros(D, D)

        # [a_i^†, a_j^†] == 0
        @test commutator(b_i |> transpose, b_j |> transpose) ≈ zeros(D, D)

        # [a_i, a_j^†] == δ_ij
        #@test commutator(b_i, b_j |> transpose) ≈ δ(i, j, D)
    end

end
