
@testset "Lattice (1D Bose-Hubbard)" begin
    for T ∈ (Float16, Float32, Float64)
        M = 3
        J = T(1)
        U = T(1/2)

        inst = Dict((i, i+1) => J for i ∈ 1:M-1)
        for i ∈ 1:M push!(inst, (i, i) => U) end

        bhg = lattice(T, inst)

        @test nv(bhg) == M
        @test ne(bhg) == M - 1

        for edge ∈ edges(bhg)
            Jg = get_prop(bhg, edge, :J)
            @test typeof(Jg) == T
            @test Jg == J
        end

        for v ∈ 1:nv(bhg)
            Ug = get_prop(bhg, v, :U)
            @test typeof(Ug) == T
            @test Ug == U
        end

        bhg_1d = chain(M, J, U, Val(:OBC))

        @test edges(bhg) == edges(bhg_1d)
        @test nv(bhg) == nv(bhg_1d)
    end
end
