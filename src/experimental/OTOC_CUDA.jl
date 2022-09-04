
"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Complex, H, ket, ϵ::Real=1E-6)
    A(du, u, p, t) = mul!(du, c .* H, u)
    ODE = ODEProblem(A, ket, (zero(eltype(τ)), τ))
    alg = RK4()
    sol = solve(ODE, alg)
    Uket = sol[end]
    Uket ./ norm(Uket)

    #=
    if τ ≈ zero(eltype(τ)) return ket end
    T = eltype(ket)

    δt = (ϵ / τ) ^ (1 / 3)
    times = fill(δt, floor(Int, τ / δt))
    tot = sum(times)
    if !(tot ≈ τ) times = [times..., abs(tot - τ)] end

    Uket = copy(ket)
    A = c .* H

    for δt ∈ times
        k1 = δt .* A * Uket
        k2 = δt .* A * (Uket .+ T(1/2) .* k1)
        k3 = δt .* A * (Uket .+ T(1/2) .* k2)
        k4 = δt .* A * (Uket .+ k3)
        Uket .+= (k1 .+ T(2) .* (k2 .+ k3) .+ k4) ./ T(6)
        Uket ./= norm(Uket)
    end
    Uket
    =#
end

function OTOC(
    times::Vector{S}, ham::BoseHubbard, i::Int, j::Int, state::State{T}, dev::Symbol, ϵ::Real=1E-6
) where {S <: Real, T <: Complex}
    @assert dev ∈ (:CPU, :GPU)

    ϕ₀ = dense(state, ham.basis)
    ai = annihilation(T, ham.basis, i)
    aj = annihilation(T, ham.basis, j)
    H = ham.H

    if dev == :GPU
        H, ai, aj = CuSparseMatrixCSC.((H, ai, aj))
        ϕ₀ = CUDA.CuArray(ϕ₀)
    end
    ψ₀ = ai * ϕ₀

    otoc = zeros(T, length(times))
    for (i, τ) ∈ enumerate(times)
        # 1. ψ := V * aj * U * ai * ϕ₀
        ψ = expv(τ, -1im, H, ψ₀, ϵ)
        ψ = expv(τ, 1im, H, aj * ψ, ϵ)

        # 2. ϕ := ai * V * aj * U * ϕ₀
        ϕ = expv(τ, -1im, H, ϕ₀, ϵ)
        ϕ = ai * expv(τ, 1im, H, aj * ϕ, ϵ)

        # 3. OTOC: <ϕ|ψ>
        otoc[i] = dot(ϕ, ψ)
    end
    otoc
end
