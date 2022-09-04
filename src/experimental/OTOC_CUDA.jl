
"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Complex, H, ket)
    A(du, u, p, t) = mul!(du, c .* H, u)
    ODE = ODEProblem(A, ket, (zero(eltype(τ)), τ))
    alg = RK4()
    sol = solve(ODE, alg, save_everystep=false)
    Uket = sol[end]
    Uket ./ norm(Uket)
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
