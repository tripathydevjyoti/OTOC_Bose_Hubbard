
@testset "OTOC" begin
    K = 4
    M, N = 3, 3

    T = Float64
    J, U = T(4/10), zero(T)

    graph = star_digraph(M)
    B = NBasis([N, N-1, N-2], M)
    H = BoseHubbard(B, J, U, graph)

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

    state = State(rand(T, K), H.basis.eig_vecs[1:K])
    otoc = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end
end
