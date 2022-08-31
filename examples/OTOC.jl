
using OTOC_Bose_Hubbard
using LightGraphs
using Plots

function bench(N, M, graph)
    K = 4

    T = Float64
    J, U = T(4/10), zero(T)

    B = NBasis.([N, N-1, N-2], M)
    H = BoseHubbard.(B, J, U, Ref(graph))

    times = [zero(T) + T(1/10) * i for i ∈ 1:100]

    state = State(rand(T, K), H[1].basis.eig_vecs[1:K])
    otoc = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end
    times, otoc
end

M = 8
N = 8
graph = star_digraph(M);

times, otoc = bench(N, M, graph);

p = plot(times, real.(otoc), title="OTOC N=$N M=$M")
#savefig(p, "./otoc.pdf")
