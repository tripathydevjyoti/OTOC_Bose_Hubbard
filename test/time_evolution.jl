using KrylovKit
using LinearAlgebra

M = N = 3
D = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))

B = Basis(M, N; constraint=:conserved_particles)

@testset "Total number of particles is concerved with time" begin
    T = Float64
    J = T(1)
    U = T(1/2)
    H = hamiltonian(M, N, J, U, :OBC)

    ket = dense_eigen_vec(B, [1, 2, 0])
    num_total = occupation(T, B)

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]
    avg = Complex{T}[]
    converged = Int[]
    
    for t ∈ times
        Uket, info = exponentiate(H, -1im * t, ket, ishermitian=true)
        push!(converged, info.converged)
        push!(avg, dot(Uket, num_total * Uket))
    end

    @test all(converged .== 1)
    @test isapprox(imag.(avg), zeros(T, length(times)), atol = 1E-14)
    @test all(T(N) .≈ real.(avg))
end
