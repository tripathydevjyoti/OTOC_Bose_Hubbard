
using OTOC_Bose_Hubbard
using LightGraphs
using Plots

function bench(N, M)
    T = Float64
    J, U = T(4/10), zero(T)

    H = BoseHubbard.([N, N-1, N-2], M, J, U, :OBC)

    times = [zero(T) + T(1/10) * i for i ∈ 1:500]

    state = State([one(T)], [fill(1, M)])
    otoc = []
    i, j = 1, 2
    for t ∈ times
        push!(otoc, OTOC(H, i, j, state, t))
    end
    times, otoc
end

M = 8
N = 8

times, otoc = bench(N, M);

p = plot(times, abs.(otoc), title="OTOC N=$N M=$M")
savefig(p, "./examples/otoc.pdf")
