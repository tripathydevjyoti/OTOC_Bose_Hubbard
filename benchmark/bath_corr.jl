include("../src_compare/OTOC_Bose_Hubbard.jl")
using .OTOC_Bose_Hubbard

using Test
using LightGraphs
using LabelledGraphs

using Plots
using PyCall

using KrylovKit
using LinearAlgebra

function expv(τ::Number, ham::BoseHubbard{T}, v::State; kwargs=()) where T
    U_dket, info = exponentiate(
        ham.H, τ, dense(v, ham.basis), ishermitian=true, tol=1E-8
    )
    @assert info.converged == 1
    U_dket
end




function bath_bartek(dim::Dims, time1::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time1)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0,1.2,50)
    state = State([one(T)], [fill(1, M)])
    J .* times, bath.(times, Ref(H), site1, site2, Ref(state))
end  

function bath_bartek2(dim::Dims, time1::Real, num_points::Int, site1::Int, site2::Int)
    T = eltype(time1)
    J, U = T(4), T(16)
    graph = hexagonal_graph(dim, J::T, U::T, :OBC)
    #graph = hex_graph(T(4),T(16))
    M = nv(graph)
    N = Int(M / 1)
    H = BoseHubbard.(NBasis.([N+1, N, N-1,N-2], M), Ref(graph))
    #times = zero(T) .+ T(time / num_points) .* collect(1:num_points)
    times = range(0,1.2,50)
    state = State([one(T)], [fill(1, M)])
    
    J .* times, bath2.(times, Ref(H), site1, site2, Ref(state))
end   



dim = (1, 1)
time1 =0.25
num_points = 40

@time Jtimes, bath_corr11 = bath_bartek2(dim,time1,num_points,4,4)
@time Jtimes, bath_corr21 = bath_bartek(dim, time1, num_points,1,4)
@time Jtimes, bath_corr31 = bath_bartek(dim, time1, num_points,4,3)
@time Jtimes, bath_corr41 = bath_bartek(dim, time1, num_points,1,3)
@time Jtimes, bath_corr51 = bath_bartek(dim, time1, num_points,4,8)
@time Jtimes, bath_corr61 = bath_bartek(dim, time1, num_points,1,8)

@time Jtimes, bath_corr71 = bath_bartek(dim, time1, num_points,3,8)
@time Jtimes, bath_corr81 = bath_bartek(dim, time1, num_points,8,3)
@time Jtimes, bath_corr91 = bath_bartek(dim, time1, num_points,8,4)

@time Jtimes, bath_corr1 = bath_bartek(dim, time, num_points,3,4)
@time Jtimes, bath_corr2 = bath_bartek(dim, time, num_points,3,8)
@time Jtimes, bath_corr3 = bath_bartek(dim, time, num_points,1,4)
@time Jtimes, bath_corr4 = bath_bartek(dim, time, num_points,4,1)

@time Jtimes, bath_corr5 = bath_bartek(dim, time, num_points,4,4)
@time Jtimes, bath_corr6 = bath_bartek(dim, time, num_points,8,8)
@time Jtimes, bath_corr7 = bath_bartek(dim, time, num_points,3,3)
@time Jtimes, bath_corr8 = bath_bartek(dim, time, num_points,4,8)
@time Jtimes, bath_corr9 = bath_bartek(dim, time, num_points,4,8)

@time Jtimes, bath_corr10 = bath_bartek(dim, time, num_points,3,1)
@time Jtimes, bath_corr11 = bath_bartek(dim, time, num_points,4,4)
@time Jtimes, bath_corr12 = bath_bartek(dim, time, num_points,8,4)

@time Jtimes, bath_corr13 = bath_bartek(dim, time, num_points,1,3)
@time Jtimes, bath_corr14 = bath_bartek(dim, time, num_points,4,4)
@time Jtimes, bath_corr15 = bath_bartek(dim, time, num_points,4,8)

@time Jtimes, bath_corr16 = bath_bartek(dim, time, num_points,4,6)
@time Jtimes, bath_corr17 = bath_bartek(dim, time, num_points,6,6)
@time Jtimes, bath_corr18 = bath_bartek(dim, time, num_points,6,8)



plot(Jtimes,[real.(bath_corr11),real.(bath_corr21),real.(bath_corr31),real.(bath_corr41)])

bath_corr=[]
for i in 1:length(Jtimes)
    push!(bath_corr,(bath_corr11[i]+ bath_corr21[i]+ bath_corr31[i]+ bath_corr41[i]+ bath_corr51[i]+ bath_corr61[i])) 
end
np = pyimport("numpy")
np.save("bath_corr_2hex_flake_ba+",bath_corr)

p = plot(Jtimes, real.(bath_corr))
xlabel!("Jt")
ylabel!(" ⟨ B(3,2)(t) B(5,1)⟩  for 2 hex flake")
savefig(p," ⟨ B(3,2)(t) B(5,1)⟩_2hex_flake.pdf" )
np = pyimport("numpy")
np.save("bath_corr_1hex_4_2",bath_corr)