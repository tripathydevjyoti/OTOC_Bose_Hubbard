using OTOC_Bose_Hubbard
using LinearAlgebra
using Plots
using LaTeXStrings
using LightGraphs
using FFTW


# compute Loschmidt Echo and Work distribution for the BH model with unit filling

function LE_1D(t::Float64)

    T = eltype(t)
    J, U = T(1), T(10)
    
    #1D parameters
    M,N=4,4
    graph = chain(M, J, U, :PBC)

    hamiltonian=BoseHubbard(N, M, J, U, graph) #BH hamiltonian
    state_i=State([one(T)], [fill(1, M)]) #initial state
    state_f=expv(-1im*t,hamiltonian,state_i) #final state
    ket=dense(state_i,NBasis(N,M))

    le=conj(dot(ket,state_f)) #complex conjugate of O(t)
    le
end

function LE_2D(t::Float64)

    T = eltype(t)
    J, U = T(1), T(10)

    #2D parameters
    dim=(1,2)
    graph= hexagonal_graph(dim, J, U, :OBC)
    M = nv(graph)
    N = Int(M / 1)

    hamiltonian=BoseHubbard(N, M, J, U, graph) #BH hamiltonian
    state_i=State([one(T)], [fill(1, M)]) #initial state
    state_f=expv(-1im*t,hamiltonian,state_i) #final state
    ket=dense(state_i,NBasis(N,M))

    le=abs(dot(ket,state_f))^2 #loschmidt echo
    le
end

function work_dis(t::Float64)

    T = eltype(t)
    J, U = T(1), T(10)
    
    #1D parameters
    M,N=4,4
    graph = chain(M, J, U, :PBC)

    hamiltonian=BoseHubbard(N, M, J, U, graph) #BH hamiltonian
    state_i=State([one(T)], [fill(1, M)]) #initial state
    state_f=expv(-1im*t,hamiltonian,state_i) #final state
    ket=dense(state_i,NBasis(N,M))

    o=conj(dot(ket,state_f)) #complex conjugate of O(t)
    o
end


interval = range(0, 1000, length=1000)
y = work_dis.(interval)
F = ifft(y)
freqs = fftfreq(length(interval))

p=plot(freqs,abs.(F),label=L"U=10")
xlabel!(L"W")
ylabel!(L"\mathcal{P}(W)")
#savefig(p, "./examples/work_pdf_1D.pdf")
p
