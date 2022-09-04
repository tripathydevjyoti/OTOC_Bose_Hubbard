export
    OTOC_ODE_CUDA

"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Number, H, ket)
    T = eltype(ket)

    num_points = 500
    Uket = copy(ket)
    dt = c * τ / num_points

    A = dt .* H
    for _ ∈ num_points
        k1 = A * Uket
        k2 = A * (Uket .+ T(1/2) .* k1)
        k3 = A * (Uket .+ T(1/2) .* k2)
        k4 = A * (Uket .+ k3)
        Uket += sum(T[1/6, 1/3, 1/3, 1/6] .* [k1, k2, k3, k4])
    end
    Uket ./ norm(Uket)
end

#=
"""
$(TYPEDSIGNATURES)
"""
function OTOC_ODE_CUDA(
    time::Vector{T}, H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State
) where {S, T <: Number}
    otoc = Complex{T}[]

    @time begin
        HN = CUSPARSE.CuSparseMatrixCSC(H[1].H)
        HN1 = CUSPARSE.CuSparseMatrixCSC(H[2].H)
        HN2 = CUSPARSE.CuSparseMatrixCSC(H[3].H)
    end

    ai_ket = dense(destroy(state, i), H[2].basis)
    ket = dense(state, H[1].basis)

    for τ ∈ time
        # 1. |x> := V * a_j * U * a_i * ket
        U_ai_ket = expv(τ, -1im, HN1, ai_ket) |> Array

        aj_U_ai_ket = dense(destroy(State(U_ai_ket, H[2].basis), j), H[3].basis)
        V_aj_U_ai_ket = expv(τ, 1im, HN2, aj_U_ai_ket) |> Array

        # 2. |y> := a_i * V * a_j * U * ket
        U_ket = expv(τ, -1im, HN, ket) |> Array
        aj_U_ket = dense(destroy(State(U_ket, H[1].basis), j), H[2].basis)

        V_aj_U_ket = expv(τ, 1im, HN1, aj_U_ket) |> Array
        ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

        # 3. OTOC: <y|x>
        push!(otoc, dot(ai_V_aj_U_ket, V_aj_U_ai_ket))
    end
    otoc
end
=#

function OTOC_ODE_CUDA(
    times::Vector{T}, ham::BoseHubbard{S}, i::Int, j::Int, state::State
) where {S, T <: Number}
    C = Complex{T}

    H = CUSPARSE.CuSparseMatrixCSC(ham.H)
    ai = CUSPARSE.CuSparseMatrixCSC(annihilation(C, ham.basis, i))
    aj = CUSPARSE.CuSparseMatrixCSC(annihilation(C, ham.basis, j))

    ket = CUDA.CuArray(dense(state, ham.basis))
    ai_ket = ai * ket

    otoc = zeros(C, length(times))
    for (i, τ) ∈ enumerate(times)
        # 1. |x> := V * a_j * U * a_i * ket
        U_ai_ket = expv(τ, -1im, H, ai_ket)
        V_aj_U_ai_ket = expv(τ, 1im, H, aj * U_ai_ket)

        # 2. |y> := a_i * V * a_j * U * ket
        U_ket = expv(τ, -1im, H, ket)
        ai_V_aj_U_ket = ai * expv(τ, 1im, H, aj * U_ket)

        # 3. OTOC: <y|x>
        otoc[i] = dot(ai_V_aj_U_ket, V_aj_U_ai_ket)
    end
    otoc
end
