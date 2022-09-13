using Test
using JeszenszkiBasis

@testset "Check Super destroy,create" begin
    basis1 = Szbasis(3,2)
    basis2 = Szbasis(3,1)
    Ψ = [1/5+0im , 0+0im ,1/5+0im, 0+0im , 0+0im, 0+0im]
    ϕ = [1/5+0im, 0+0im, 1/5+0im]
    site = 1
    @test super_destroy(Ψ, basis1, basis2, site) == [0.2,0,0]
    @test super_create(ϕ, basis2, basis1, site) == [0.2,0,0,0.2,0,0]

end


