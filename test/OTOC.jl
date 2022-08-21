
@testset "OTOC" begin
    M, N = 3, 3

    ϵ = 1E-10
    T = Float64
    J = T(4/10)
    U = zero(T)

    H = BoseHubbard(N, M, J, U, :OBC).H

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

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
