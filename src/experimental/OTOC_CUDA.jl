export
    OTOC_CUDA

"""
$(TYPEDSIGNATURES)
"""
function expv(τ::Real, c::Number, ham::BoseHubbard{T}, state::State) where T
    f(du, u, p, t) = mul!(du, c .* ham.H, u)
    ket = convert.(Complex{T}, dense(state, ham.basis))
    sol = solve(ODEProblem(f, ket, τ), Tsit5())
    sol(τ)
end

function OTOC_CUDA(
    H::Vector{BoseHubbard{S}}, i::Int, j::Int, state::State, time::Vector{T}
) where {S, T <: Number}
    otoc = Complex{T}[]

    HN = CuSparseMatrixCSC(H[1].H)
    HN1 = CuSparseMatrixCSC(H[2].H)
    HN2 = CuSparseMatrixCSC(H[3].H)

    ai_ket = CUDA.CuArray(dense(destroy(state, i), H[1].basis))
    ket = CUDA.CuArray(dense(state, i))

    for t ∈ time
        # 1. compute |x> := V * a_j * U * a_i * ket
        U_ai_ket = expv(t, HN1, ai_ket) |> Array
        x = CUDA.CuArray(destroy(State(U_ai_ket, H[2].basis), j))
        V_aj_U_ai_ket = expv(t, HN2, x) |> Array

        # 2. compute |y> := a_i * V * a_j * U * ket
        U_ket = expv(t, HN, ket) |> Array
        V_aj_U_ket = expv(τ, HN1, destroy(State(U_ket, H[1].basis), j))
        ai_V_aj_U_ket = dense(destroy(State(V_aj_U_ket, H[2].basis), i), H[3].basis)

        # 3. compute OTOC: <y|x>
        push!(otoc, dot(ai_V_aj_U_ket, V_aj_U_ai_ket))
    end
    otoc
end
