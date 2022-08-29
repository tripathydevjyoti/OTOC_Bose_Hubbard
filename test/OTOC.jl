
@testset "OTOC" begin
    num_points = 100
    K = 4
    M, N = 3, 3

    T = Float64
    J, U = T(4/10), zero(T)

    graph = star_digraph(M)

    B = NBasis([N, N-1, N-2], M)
    H = BoseHubbard(B, J, U, graph)

    B2 = NBasis.([N, N-1, N-2], Ref(M))
    H2 = BoseHubbard.(B2, Ref(J), Ref(U), Ref(graph))

    times = [zero(T) + T(1/10) * i for i ∈ 1:num_points]

    state = State(rand(T, K), H2[1].basis.eig_vecs[1:K])
    state_2 = State(rand(T, K), H2[1].basis.eig_vecs[1:K])

    otoc = []
    otoc_2 = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
        push!(otoc_2, OTOC(H2, i, j, state, t))
    end

    @test length(otoc) == length(otoc_2) == num_points
    @test otoc ≈ otoc_2
end
