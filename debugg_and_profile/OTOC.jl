
using OTOC_Bose_Hubbard
using LightGraphs
using Profile, PProf
using FlameGraphs

function bench(N, M, K, graph)
    T = Float64
    B = NBasis.([N, N-1, N-2], Ref(M))
    OTOC(
        BoseHubbard.(B, T(4/10), zero(T), Ref(graph)),
        1, 2,
        State(rand(T, K), B[1].eig_vecs[1:K]),
        one(T)
    )
    0
end

M = 10
N = 10
K = 10
graph = star_digraph(M)

bench(N, M, K, graph)
#@profile bench(N, M, K, graph)
@time Profile.Allocs.@profile sample_rate = 1 bench(N, M, K, graph)
PProf.Allocs.pprof()

#pprof(flamegraph(); webhost = "localhost", webport = 51456)
