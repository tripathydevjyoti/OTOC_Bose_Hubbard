
"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Complex, H, ket)
    z = zero(eltype(τ))
    if τ ≈ z return ket end

    A(du, u, p, t) = mul!(du, c .* H, u)
    ODE = ODEProblem(A, ket, (z, τ))
    alg = RK4()
    #alg = ETDRK4(krylov=true, m=30)

    sol = solve(
        ODE, alg, save_everystep=false, reltol=1e-8, abstol=1e-8
    )

    Uket = sol[end]
    Uket ./ norm(Uket)

#=
    T = eltype(ket)
    δt = (1e-6 / τ) ^ (1 / 3)
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
    times::Vector{S}, ham::BoseHubbard, i::Int, j::Int, state::State{T}, dev::Symbol
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
        ψ = expv(τ, -1im, H, ψ₀)
        ψ = expv(τ, 1im, H, aj * ψ)

        # 2. ϕ := ai * V * aj * U * ϕ₀
        ϕ = expv(τ, -1im, H, ϕ₀)
        ϕ = ai * expv(τ, 1im, H, aj * ϕ)

        # 3. OTOC: <ϕ|ψ>
        otoc[i] = dot(ϕ, ψ)
    end
    otoc
end
