export
    OTOC_bose_bubbard

"""
$(TYPEDSIGNATURES)
"""
function OTOC_bose_bubbard(
    H::S, a::S, b:::S,
    ket::Array{Complex{T}},
    time::T;
    kwargs
) where {T <: Real, S <: SparseArrays}

    U, infoU = exponentiate(H, -1im * time, ket; (kwargs..., ishermitian=true)...)
    V, infoV = exponentiate(H,  1im * time, ket; (kwargs..., ishermitian=true)...)

    dot(a * V * b * U * ket, V * b * U * a * ket)
end
