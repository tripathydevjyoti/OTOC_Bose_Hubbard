dim(N::Int, M::Int) = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))

@testset "Dense Hamiltonian" begin
    M = N = 3
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

    T = Float64
    J = T(1)
    U = T(1/2)
    H = hamiltonian(N, M, J, U, :OBC)

    @test size(H) == (D, D)
    @test transpose(H) == H
end

@testset "Toy model with 2 sites and 1 particle" begin
    M, N = 2, 1
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

    ϵ = 1E-10
    T = Float64
    J = T(4/10)
    U = zero(T)

    H = hamiltonian(N, M, J, U, :OBC)

    @test size(H) == (D, D)
    @test transpose(H) == H

    evals, eigen_vec = eigen(Hermitian(H |> Array))

    @test evals ≈ [-J, J]

    w0 = dense_eigen_vec(B, [1, 0])
    w1 = dense_eigen_vec(B, [0, 1])

    v0 = (w0 .+ w1) ./ sqrt(2)
    v1 = (w0 .- w1) ./ sqrt(2)

    @test eachcol(eigen_vec) |> collect ≈ [-v0, -v1]
end
