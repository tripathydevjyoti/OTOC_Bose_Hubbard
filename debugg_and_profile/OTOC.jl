
using OTOC_Bose_Hubbard
using LightGraphs
using Profile, PProf
using FlameGraphs

function bench(N, M, K, graph)
    T = Float64
    B = NBasis([N, N-1, N-2], M)
    println("Dim: ", B.dim)
    OTOC(
        BoseHubbard(B, T(4/10), zero(T), graph),
        1, 2,
        State(rand(T, K), B.eig_vecs[1:K]),
        one(T)
    )
    0
end

M = 8
N = 8
K = 10
graph = star_digraph(M)

bench(N, M, K, graph)
@profile bench(N, M, K, graph)

pprof(flamegraph(); webhost = "localhost", webport = 54326)
