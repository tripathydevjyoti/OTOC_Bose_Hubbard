
@testset "OTOC" begin
    num_points = 100
    M, N = 3, 3

    T = Float64
    J, U = T(4/10), zero(T)

    graph = star_digraph(M)

    B = NBasis([N, N-1, N-2], M)
    H = BoseHubbard(B, J, U, graph)

    B2 = NBasis.([N, N-1, N-2], M)
    H2 = BoseHubbard.(B2, J, U, Ref(graph))

    B3 = Basis(N, M)
    H3 = BoseHubbard(B3, J, U, graph)

    times = [zero(T) + T(1/10) * (i-1) for i ∈ 1:num_points]

    i, j = 1, 2
    eigs = [[1, 1, 1], [1, 2, 0], [2, 1, 0]]
    coeff = rand(T, length(eigs))
    coeff ./= sqrt(sum(abs.(coeff) .^ 2))

    state = State(coeff, eigs)

    otoc = OTOC.(times, Ref(H), i, j, Ref(state))
    otoc_2 = OTOC.(times, Ref(H2), i, j, Ref(state))
    otoc_3 = OTOC.(times, Ref(H3), i, j, Ref(state))

    ni_nj = getindex.(state.eig_vecs, i) .* getindex.(state.eig_vecs, j)

    @test all(sum.(eigs) .== N)
    @test length(otoc) == length(otoc_2) == num_points
    @test otoc ≈ otoc_2 ≈ otoc_3
    @test times[1] ≈ zero(T)
    @test real(otoc[1]) ≈ sum(abs.(state.coeff) .^ 2 .* ni_nj)
    @test imag(otoc[1]) ≈ zero(T)
end
