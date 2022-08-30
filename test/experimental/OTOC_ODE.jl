
@testset "OTOC via ODE" begin
    num_points = 100
    K = 4
    M, N = 3, 3

    T = Float64
    J, U = T(4/10), zero(T)

    graph = star_digraph(M)

    B = NBasis.([N, N-1, N-2], M)
    H = BoseHubbard.(B, J, U, Ref(graph))

    times = [zero(T) + T(1/10) * i for i ∈ 1:num_points]

    state = State(rand(T, K), H[1].basis.eig_vecs[1:K])

    otoc = Complex{T}[]

    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end
    otoc_ode = OTOC_ODE(H, i, j, state, times)

    @test length(otoc) == length(otoc_ode) == num_points
    @test isapprox(otoc, otoc_ode, atol = 1E-2)
end
