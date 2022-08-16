M = N = 3

B_ref = [
    [3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0],
    [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1],
    [0, 1, 2], [0, 0, 3]
]

@testset "Basis" begin
    B = Basis(M, N)
    @test eltype(B) == Uint
    @test length.(B) == M
    @test Set(B) == Set(B_ref)
    @test sum.(B) == fill(N, length(B))
    println(B)
end
