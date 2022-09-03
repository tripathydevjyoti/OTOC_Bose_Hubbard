
using OTOC_Bose_Hubbard
using LightGraphs
using CUDA
using CUDA.CUSPARSE
using KrylovKit

function bench(N, M, graph)
    T = Float64
    S = Complex{T}

    B = NBasis([N, N-1, N-2], M)
    ham = BoseHubbard(B, S(4/10), zero(S), graph)

    v = CUDA.rand(S, ham.basis.dim)
    @time H = CuSparseMatrixCSC(ham.H)
    println(typeof(H))

    @time Ux, info = exponentiate(H, -1im, v, ishermitian=true)
    println(typeof(Ux))
    @assert info.converged == 1
end

M = 5
N = 5
graph = star_digraph(M);

bench(N, M, graph);
bench(N, M, graph);

nothing
