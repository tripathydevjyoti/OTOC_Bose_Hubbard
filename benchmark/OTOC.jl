
using OTOC_Bose_Hubbard
using LightGraphs

function bench(N, M, graph)
    K = 4

    T = Float64
    J, U = T(4/10), zero(T)

    B = NBasis.([N, N-1, N-2], Ref(M))
    H = BoseHubbard.(B, Ref(J), Ref(U), Ref(graph))

    #println("Dim: ", B.dim)

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

    state = State(rand(T, K), H[1].basis.eig_vecs[1:K])
    otoc = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end
    otoc
end

M = 10
N = 10
graph = star_digraph(M);

@time otoc = bench(N, M, graph);
@time otoc = bench(N, M, graph);

nothing
