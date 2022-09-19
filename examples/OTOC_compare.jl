
using OTOC_Bose_Hubbard
using LightGraphs
using Plots

#=
function bench(dim::Dims, time::Real, ij::Dims, num_points::Int)
    T = eltype(time)
    J, U = T(4), T(16)

    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard(NBasis([N, N-1, N-2], M), graph)

    times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    state = State([one(Complex{T})], [fill(1, M)])

    J .* times, OTOC(times, H, ij..., state, :GPU)
end
=#

function bench(dim::Dims, time::Real, ij::Dims, num_points::Int)
    T = eltype(time)
    J, U = T(4), T(16)

    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N, N-1, N-2], M), Ref(graph))

    times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    state = State([one(T)], [fill(1, M)])

    J .* times, OTOC.(times, Ref(H), ij..., Ref(state))
end

dim = (1, 1)
time = 1.0
num_points = 50

ij = (6, 1)
@time Jtimes, otoc = bench(dim, time, ij, num_points)
p = plot(Jtimes, abs.(otoc), title="OTOC for hex $dim", label="$ij")

ij = (2, 1)
@time Jtimes, otoc = bench(dim, time, ij, num_points)
plot!(p, Jtimes, abs.(otoc), label="$ij")

savefig(p, "./examples/otoc_compare.pdf")
