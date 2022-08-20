
using OTOC_Bose_Hubbard

function bench(M::Int, N::Int)
    J, U = 1.0, 0.5
    BoseHubbard(N, M, J, U, :OBC)
end

M = N = 13
@time bench(M, N)
@time bench(M, N)
nothing
