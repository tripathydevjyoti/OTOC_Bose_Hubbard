
using OTOC_Bose_Hubbard
using LightGraphs
using Plots

function bench(::Type{T}, dim::Dims, time::Real, num_points::Int) where T
    J, U = T(4), T(16)

    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard(NBasis([N, N-1, N-2], M), graph)

    times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    state = State([one(Complex{T})], [fill(1, M)])

    times, OTOC_ODE_CUDA(times, H, 1, 2, state)
end

dim = (1, 2)
time = 5.0
num_points = 200
@time times, otoc = bench(Float64, dim, time, num_points)

p = plot(times, abs.(otoc), title="OTOC for hex $dim")
savefig(p, "./examples/otoc_cuda_hex$dim.pdf")
