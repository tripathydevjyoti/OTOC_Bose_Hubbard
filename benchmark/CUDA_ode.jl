
using OTOC_Bose_Hubbard
using LightGraphs
using LinearAlgebra
using CUDA
using CUDA.CUSPARSE
using DifferentialEquations

function bench(N, M, graph)
    T = Float64
    τ = one(T)

    @time begin
        B = NBasis(N, M)
        ham = BoseHubbard(B, T(4/10), zero(T), graph)
        H = CuSparseMatrixCSC(ham.H)
    end
    v = CUDA.rand(T, ham.basis.dim)

    f(du, u, p, t) = mul!(du, H, u)
    prob = ODEProblem(f, v,  τ)
    @time sol = solve(prob, Tsit5())
    println(typeof(sol(τ)))
end

M = 10
N = 10
graph = star_digraph(M);

bench(N, M, graph);
bench(N, M, graph);

nothing
