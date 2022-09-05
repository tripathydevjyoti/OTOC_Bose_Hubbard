
@testset "OTOC via ODE" begin
    ϵ = 1e-6
    T = Float64

    dim = (1, 1)
    time = T(1)
    num_points = 2

    J, U = one(T), zero(T)

    M = 3
    graph = chain(M, J, U, :OBC)
    #graph = hexagonal_graph(dim, J, U, :OBC)
    #M = nv(graph)
    N = M #Int(M / 1)

    B = NBasis.([N, N-1, N-2], M)
    H = BoseHubbard.(B, J, U, Ref(graph))

    times = zero(T) .+ T(time / num_points) .* collect(0:num_points-1)
    state = State([one(Complex{T})], [fill(1, M)])

    i, j = 1, 2
    @time otoc = OTOC.(times, Ref(H), i, j, Ref(state))
    @time otoc_ode = OTOC_ODE(times, H, i, j, state)

    @test length(otoc) == length(otoc_ode) == num_points
    @test isapprox(otoc, otoc_ode, atol = ϵ)

    if J ≈ U ≈ zero(T) all(otoc .≈ zero(Complex{T})) end

    if J ≈ zero(T)
        otoc_diag = Complex{T}[]

        ket = dense(state, H[1].basis)
        ai_ket = destroy(state, i)

        for τ ∈ times
            U_ai_ket = exp(-1im * τ .* Array(H[2].H)) * dense(ai_ket, H[2].basis)
            x = destroy(State(U_ai_ket, H[2].basis), j)
            a = exp(1im * τ .* Array(H[3].H)) * dense(x, H[3].basis)

            U_ket = exp(-1im * τ .* Array(H[1].H)) * ket
            y = destroy(State(U_ket, H[1].basis), j)
            V_aj_U_ket = exp(1im * τ .* Array(H[2].H)) * dense(y, H[2].basis)
            b = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

            push!(otoc_diag, dot(b, a))
        end
        @test otoc ≈ otoc_diag
    end
end
