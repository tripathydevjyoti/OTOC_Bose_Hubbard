using Test
using JeszenszkiBasis



function superposition(
    Ψ ::Vector{Complex{T}},
    basis ::AbstractSzbasis
)   where T<: Real
    Ψₙ = []
    β = Complex[]

    Cₙ =broadcast(abs, Ψ)
    
    for(i,x) in enumerate(Cₙ)
        
        
        if x > 0.001
            push!(Ψₙ,basis[i])
            push!(β,Ψ[i])
        end
        
    end

        


    
    return β,Ψₙ 
    
end



@testset "Checking if superposition is working fine" begin
    basis = Szbasis(3,1) # 3 site 1 boson model 
    Ψ = [1/5 ,4/5 ,0]
    Ψ = complex(Ψ)
    
    coeff_list,states_list = superposition(Ψ , basis)
    #coeff_list = broadcast(abs , coeff_list)
    @test broadcast(abs, coeff_list) == [1/5,4/5]
    @test states_list == [[1,0,0], [0,1,0]]
end
