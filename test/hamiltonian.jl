@testset "Hamiltonian" begin
    M = N = 3
    B = Basis(N, M)

    for T ∈ (Float16, Float32, Float64)
        J, U = T(1), T(0)
        for bndr ∈ (:OBC, :PBC)
            ham = BoseHubbard{T}(B, chain(M, J, U, bndr)) #BoseHubbard(N, M, J, U, bndr)
            #@test size(ham.H) == (B.dim, B.dim)
            @test conj(transpose(ham.H)) == ham.H
            @test ham.H ≈ BoseHubbard{T}(B, chain(M, J, U, bndr), 0).H
        end
    end
end

@testset "Toy model with 2 sites and 1 particle" begin
    M, N = 2, 1

    B = NBasis(N, M)

    ϵ = 1E-10
    T = Float64
    J, U = T(4/10), zero(T)

    H = BoseHubbard(N, M, J, U, :OBC).H

    @test size(H) == (B.dim, B.dim)
    @test conj(transpose(H)) == H

    evals, eigen_vec = eigen(Hermitian(Array(H)))

    @test evals ≈ [-J, J]

    w0 = dense([1, 0], B)
    w1 = dense([0, 1], B)

    v0 = (w0 .+ w1) ./ sqrt(2)
    v1 = (w0 .- w1) ./ sqrt(2)

    @test eachcol(eigen_vec) |> collect ≈ [-v0, -v1]
end
