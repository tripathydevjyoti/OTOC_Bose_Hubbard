@testset "Operate (destroy and create)" begin
    @test operate([3, 0, 0], 1, :destroy) == [2, 0 ,0]
    @test operate([3, 0, 0], 2, :create) == [3, 1 ,0]
    @test destroy_and_create([2, 0, 1], 1, 2) == [1, 1, 1]
end

@testset "Operate (destroy and create)" begin
    M = N = 3
    D = (M + 1) ^ N
    B = Basis(M, N, constraint=:none)

    for T ∈ (Float16, Float32, Float64)
        for i ∈ 1:M, j ∈ 1:M
            a_i = annihilation(T, B, i) |> Array
            ap_j = creation(T, B, i) |> Array

            @test tr(a_i) == tr(ap_j) ≈ zero(T)
            @test commutator(a_i, ap_j |> transpose) ≈ zeros(T, D, D)
            @test commutator(ap_j, a_i |> transpose) ≈ zeros(T, D, D)

            if i == j
                @test a_i == ap_j |> transpose
                @test ap_j * a_i ≈ occupation(T, B, i) |> Array
            else
                @test commutator(a_i, ap_j) ≈ zeros(T, D, D)
            end
        end
    end
end
