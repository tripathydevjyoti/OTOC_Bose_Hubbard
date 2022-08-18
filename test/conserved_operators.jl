using LinearAlgebra
using SparseArrays

M = N = 3
D = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))

B = Basis(M, N; constraint=:conserved_particles)

@testset "Sparse occupation operators" begin
    for T ∈ (Float16, Float32, Float64)
        for i ∈ 1:M
            n_i = occupation(T, B, i)
            @test eltype(n_i) == T
            @test isdiag(n_i)
            @test transpose(n_i) == n_i
            @test size(n_i) == (D, D)
        end
        @test occupation(T, B) == sum(occupation(B, i) for i ∈ 1:M)
    end
end
