
using OTOC_Bose_Hubbard
using LightGraphs

function bench(::Type{T}, dim::Dims, time::Real, num_points::Int) where T
    J, U = T(4), T(16)

    @time graph = hexagonal_graph(dim, J::T, U::T, :OBC)

    M = nv(graph)
    N = Int(M / 2)
    println(dim, " => ",  M, " ", N, " ", ne(graph))

    @time H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    state = State([one(Complex{T})], [fill(1, M)])
    times, OTOC_ODE_CUDA(times, H, 1, 2, state)
end

dim = (2, 1)
time = 2.0
num_points = 1
@time times, otoc = bench(Float64, dim, time, num_points)

nothing
