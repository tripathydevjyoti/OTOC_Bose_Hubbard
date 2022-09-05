
using OTOC_Bose_Hubbard
using LightGraphs
using Plots

function bench(dim::Dims, time::Real, num_points::Int)
    T = eltype(time)
    J, U = T(4), T(16)

    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    state = State([one(T)], [fill(1, M)])

    J .* times, OTOC.(times, Ref(H), 1, 2, Ref(state))
end

dim = (1, 2)
time = 2.0
num_points = 100
@time Jtimes, otoc = bench(dim, time, num_points)

p = plot(Jtimes, abs.(otoc), title="OTOC for hex $dim")
savefig(p, "./examples/otoc_hex$dim.pdf")
