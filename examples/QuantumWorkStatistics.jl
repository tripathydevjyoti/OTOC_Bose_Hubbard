using OTOC_Bose_Hubbard
using LinearAlgebra
using Plots
using LaTeXStrings

# compute Loschmidt Echo for 1D BH model with unit filling

M,N=4,4

function LE(t::Float64)

    T = eltype(t)
    J, U = T(1), T(1)

    graph = chain(M, J, U, :PBC)
    hamiltonian=BoseHubbard(N, M, J, U, graph)

    state_i=State([one(T)], [fill(1, M)])
    state_f=expv(-1im*t,hamiltonian,state_i)

    ket=dense(state_i,NBasis(N,M))
    le=abs(dot(ket,state_f))^2
    le
end

interval = range(0, 10, length=1000)
y = LE.(interval)

plot(interval, y)
xlabel!(L"t")
ylabel!(L"L(t)")
