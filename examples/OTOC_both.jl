
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
    H2 = BoseHubbard(NBasis([N, N-1, N-2], M), graph)

    times = zero(T) .+ T(time / num_points) .* collect(0:num_points-1)
    state = State([one(Complex{T})], [fill(1, M)])

    @time otoc = OTOC.(times, Ref(H), 1, 2, Ref(state))
    @time otoc_cuda = OTOC(times, H2, 1, 2, state, :GPU)

    J .* times, otoc, otoc_cuda
end

dim = (1, 1)
time = 5.0
num_points = 100
Jtimes, otoc, otoc_cuda = bench(dim, time, num_points)

p = plot(Jtimes, abs.(otoc), title="OTOC for hex $dim", label="|OTOC|")
plot!(p, Jtimes, abs.(otoc_cuda), label="|OTOC-CUDA|")
savefig(p, "./examples/otoc_both_hex$dim.pdf")
