
K = 6

SB = [
    [3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0],
    [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1],
    [0, 1, 2], [0, 0, 3]
]

@testset "Basis" begin
    for (M, N) ∈ ((m, n) for m ∈ 1:K, n ∈ 1:K)
        B = Basis(N, M)
        NB = NBasis(N, M)
        I = B.sub_basis_indices

        @test length(B.eig_vecs) == length(B.tags) == (N + 1) ^ M == B.dim
        @test eltype(B.eig_vecs) == Vector{Int} == eltype(NB.eig_vecs)
        @test all(length.(B.eig_vecs) .== M) && all(length.(NB.eig_vecs) .== M)
        @test length(I) == sub_basis_dim(N, M) == length(NB.eig_vecs) == NB.dim
        @test Set(B.eig_vecs[I]) == Set(NB.eig_vecs)
        if M == N == 3 @test Set(SB) == Set(NB.eig_vecs) end
        @test all(sum.(B.eig_vecs[I]) .== N) && all(sum.(NB.eig_vecs) .== N)
        @test B.N == N == NB.N && B.M == M == NB.M
        @test issorted(B.tags) && issorted(NB.tags)
        @test length(I) == NB.dim
    end
end

@testset "Superposition" begin
    M = N = 3
    NB = NBasis(N, M)

    for T ∈ (Float16, Float32, Float64)
        coeff = rand(Complex{T}, NB.dim)
        state = State(coeff, NB.eig_vecs)

        @test eltype(state) == eltype(coeff)
        @test length(state.eig_vecs) == NB.dim == length(state.coeff)
    end

    vec = sprandn(NB.dim, 0.25) |> Array
    idx = findall(!iszero, vec)
    state = State(vec[idx], NB.eig_vecs[idx])
    @test dense(state, NB) ≈ vec
end
