@testset "Hamiltonian" begin
    M = N = 3
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

    T = Float64
    J = T(1)
    U = T(1/2)
    ham = BoseHubbard(N, M, J, U, :OBC)

    #@test ham.basis == B
    #@test ham.lattice == chain(M, J, U, Val(:OBC))
    @test size(ham.H) == (D, D)
    @test transpose(ham.H) == ham.H
end

@testset "Toy model with 2 sites and 1 particle" begin
    M, N = 2, 1
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

    ϵ = 1E-10
    T = Float64
    J = T(4/10)
    U = zero(T)

    H = BoseHubbard(N, M, J, U, :OBC).H

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
