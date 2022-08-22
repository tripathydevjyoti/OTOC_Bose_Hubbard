
@testset "OTOC" begin
    num_points = 100
    K = 4
    M, N = 3, 3

    T = Float64
    J, U = T(4/10), zero(T)

    graph = star_digraph(M)
    B = NBasis([N, N-1, N-2], M)
    H = BoseHubbard(B, J, U, graph)

    times = [zero(T) + T(1/10) * i for i ∈ 1:num_points]

    state = State(rand(T, K), B.eig_vecs[1:K])
    otoc = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end

    @test length(otoc) == num_points
end
