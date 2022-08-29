
using OTOC_Bose_Hubbard
using LightGraphs
using CUDA
using CUDA.CUSPARSE
using KrylovKit

function bench(N, M, graph)
    T = Float64
    B = NBasis([N, N-1, N-2], M)
    ham = BoseHubbard(B, T(4/10), zero(T), graph)

    v = CUDA.rand(T, ham.basis.dim)
    H = CuSparseMatrixCSC(ham.H)

    @time Ux, info = exponentiate(H, -1im, v, ishermitian=true)
    println(typeof(Ux))
    @assert info.converged == 1
end

M = 8
N = 8
graph = star_digraph(M);

bench(N, M, graph);
bench(N, M, graph);

nothing
