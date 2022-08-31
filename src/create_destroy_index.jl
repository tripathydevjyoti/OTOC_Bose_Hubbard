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

export destroy
export create