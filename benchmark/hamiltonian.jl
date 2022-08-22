
using OTOC_Bose_Hubbard
using LightGraphs

function bench(N, M, graph)
    T = Float64
    @time B = NBasis([N, N-1, N-2], M)
    @time H = BoseHubbard(B, T(4/10), zero(T), graph)
    H
end

M = 10
N = 10
graph = star_digraph(M);

bench(N, M, graph);
H = bench(N, M, graph);

println(H.basis.dim)
nothing
