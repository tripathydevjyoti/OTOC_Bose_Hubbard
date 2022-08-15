using Test


function create(
    ket::Array{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] += 1
    return nket
end

function destroy(
    ket::Array{Int64,1},
    i::Int64
    )
    nket = copy(ket)
    nket[i] -= 1
    return nket
end


@testset "Checking if functions 'create' and 'destory' are working fine" begin

@test create([5,0,0,0,0,0],1) == [6,0,0,0,0,0]
    
end 