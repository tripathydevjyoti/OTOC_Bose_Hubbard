using OTOC_Bose_Hubbard
using LinearAlgebra
using Plots

M, N = 4, 4
J = 1.0

function energy_levels(U::Real)
    graph = chain(M, J, U, :PBC)
    Ham = BoseHubbard(N, M, J, U, graph).H
    eigen(Array(Ham)).values
end

function LE(t::Real)
    Ui, Uf = 0.0, 1.0
    graphi = chain(M, J, Ui,:PBC)
    graphf = chain(M, J, Uf,:PBC)
    Hi = BoseHubbard(N, M, J, Ui, graphi).H
    Hf = BoseHubbard(N, M, J, Uf, graphf).H
    ei, vi = eigen(Array(Hi))
    ef, vf = eigen(Array(Hf))
    #f = dot(vi, expv(t, Hf, vi))
    f = sum((abs(dot(vf[k],vi[1]))^2)*exp(-1im*ef[k]*t) for k ∈ 1:length(ef))
    abs(f) ^ 2
end

interval = range(0, 20, length=200)
y = LE.(interval)
#plot(interval, y)
