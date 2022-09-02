
using OTOC_Bose_Hubbard
using LightGraphs

function bench(N, M, graph, num_points=100)
    K = 4

    T = Float64
    J, U = T(4/10), zero(T)

    B = NBasis.([N, N-1, N-2], M)
    H = BoseHubbard.(B, J, U, Ref(graph))

    times = zero(T) .+ T(1/10) .* collect(1:num_points)

    state = State(rand(T, K), H[1].basis.eig_vecs[1:K])
    OTOC.(times, Ref(H), 1, 3, Ref(state))
end

M = 10
N = 10
graph = star_digraph(M);

@time otoc = bench(N, M, graph);
@time otoc = bench(N, M, graph);

nothing
