@testset "Total number of particles is concerved with time" begin
    M = N = 3
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

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

@testset "Toy model with 2 sites and 1 particle" begin
    M = 2
    N = 1
    D = dim(N, M)

    B = Basis(M, N; constraint=:conserved_particles)

    ϵ = 1E-10
    T = Float64
    J = T(4/10)
    U = zero(T)

    H = hamiltonian(N, M, J, U, :OBC)

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

    n_1 = occupation(T, B, 1)
    n_2 = occupation(T, B, 2)

    converged = Int[]
    n1, n2 = Complex{T}[], Complex{T}[]

    ket = dense_eigen_vec(B, [1, 0])
    for t ∈ times
        Uket, info = exponentiate(H, -1im * t, ket, ishermitian=true)
        push!(converged, info.converged)
        push!(n1, dot(Uket, n_1 * Uket))
        push!(n2, dot(Uket, n_2 * Uket))
    end

    @test all(converged .== 1)

    @test isapprox(imag.(n1), zeros(T, length(times)), atol = ϵ)
    @test isapprox(imag.(n2), zeros(T, length(times)), atol = ϵ)

    @test isapprox(cos.(J .* times) .^ 2, real(n1), atol = ϵ)
    @test isapprox(sin.(J .* times) .^ 2, real(n2), atol = ϵ)
end
