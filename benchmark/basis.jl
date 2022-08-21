
using OTOC_Bose_Hubbard
using Combinatorics

function bench(M::Int, N::Int)
    J, U = 1.0, 0.5
    BoseHubbard(N, M, J, U, :OBC)
end

M = N = 6
@time bench(M, N)
@time bench(M, N)

@time C2 = vec(collect.(Iterators.product(fill(collect(0:N), M)...)))

@time B1 = collect(multiexponents(M, N))
@time B2 = map(A -> [sum(A .== i) for i ∈ 1:M], with_replacement_combinations(1:M, N))
@assert B1 == B2

nothing
