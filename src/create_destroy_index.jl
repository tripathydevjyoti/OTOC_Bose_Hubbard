using JeszenszkiBasis

function create(
    ket::AbstractArray{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] += 1
    return nket
end


function destroy(
    ket::AbstractArray{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] -= 1
        
    
    
    
    return nket
end



function find_index(
    inp_state::Array{Int64,1},
    basis::AbstractSzbasis
    )
    index=0
    
    for i in 1:length(basis)
        if inp_state == basis[i]
            index =i
            break
        end

    end
    

    return index
end

function change_basis(
    ket::Vector{Int64}, basis::AbstractSzbasis
    )
    index = find_index(ket, basis)
    track = zeros(length(basis))
    track[index] = 1

    return complex(track)

end

export change_basis
export destroy
export create
export find_index