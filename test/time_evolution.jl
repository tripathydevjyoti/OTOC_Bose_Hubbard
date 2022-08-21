@testset "Total number of particles is concerved with time" begin
    M = N = 3

    for T ∈ (Float16, Float32, Float64, )
        J = T(1)
        U = T(1/2)

        B = NBasis(N, M)

        ket = dense([1, 2, 0], B)
        num_total = occupation(T, B)

        times = [zero(T) + T(1/10) * i for i ∈ 1:100]

        for graph ∈ (path_graph(M), path_digraph(M), star_digraph(M), turan_graph(M, 2))
            avg = Complex{T}[]
            converged = Int[]

            H = BoseHubbard(B, J, U, graph).H
            for t ∈ times
                Uket, info = exponentiate(H, -1im * t, ket, ishermitian=true)
                push!(converged, info.converged)
                push!(avg, dot(Uket, num_total * Uket))
            end

            @test all(converged .== 1)
            @test isapprox(imag.(avg), zeros(T, length(times)), atol = 1E-14)
            @test all(T(N) .≈ real.(avg))
        end
    end
end

@testset "Toy model with 2 sites and 1 particle" begin
    M, N = 2, 1
    B = NBasis(N, M)

    ϵ = 1E-10
    T = Float64
    J = T(4/10)
    U = zero(T)

    H = BoseHubbard(N, M, J, U, :OBC).H

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

    n_1 = occupation(T, B, 1)
    n_2 = occupation(T, B, 2)

    converged = Int[]
    n1, n2 = Complex{T}[], Complex{T}[]

    ket = dense([1, 0], B)
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
