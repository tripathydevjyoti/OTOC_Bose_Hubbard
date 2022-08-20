export
    OTOC_bose_bubbard

"""
$(TYPEDSIGNATURES)
a * V * b * U * ket, V * b * U * a * ket
"""
function OTOC_bose_bubbard(
    H::S, a::S, b:::S,
    ket::Array{Complex{T}},
    time::T;
    kwargs
) where {T <: Real, S <: SparseArrays}
    args = (kwargs..., ishermitian=true)

    Ua_ket, infoU = exponentiate(H, 1im * time, a * ket, args...)
    @assert infoU.concerved == 1
    VbUa_ket, infoV = exponentiate(H, -1im * time, b * Ua_ket, args...)
    @assert infoV.concerved == 1

    U_ket, infoU = exponentiate(H, 1im * time, ket, args...)
    @assert infoU.concerved == 1
    VbU_ket, infoV = exponentiate(H, -1im * time, b * U_ket, args...)
    @assert infoV.concerved == 1

    dot(a * VbU_ket, VbUa_ket)
end
