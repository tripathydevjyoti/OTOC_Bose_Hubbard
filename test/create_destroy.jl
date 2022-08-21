
commutator(A, B) = A * B .- A * B

@testset "Operate (destroy and create)" begin
    @test operate([3, 0, 0], 1, :destroy) == [2, 0 ,0]
    @test operate([3, 0, 0], 2, :create) == [3, 1 ,0]
    @test destroy_and_create([2, 0, 1], 1, 2) == [1, 1, 1]
end

@testset "Operate (destroy and create)" begin
    M = N = 3
    B = Basis(M, N)
    I = B.sub_basis_indices
    D = sub_basis_dim(N, M)

    for T ∈ (Float16, Float32, Float64)
        for i ∈ 1:M, j ∈ 1:M
            a_i = annihilation(T, B, i) |> Array
            ap_j = creation(T, B, i) |> Array

            @test tr(a_i) == tr(ap_j) ≈ zero(T)
            @test commutator(a_i, ap_j |> transpose) ≈ zeros(T, B.dim, B.dim)
            @test commutator(ap_j, a_i |> transpose) ≈ zeros(T, B.dim, B.dim)

            @test a_i[I, I] ≈ ap_j[I, I] ≈ zeros(T, D, D)

            if i == j
                @test a_i ≈ ap_j |> transpose
                @test ap_j * a_i ≈ occupation(T, B, i) |> Array
            else
                @test commutator(a_i, ap_j) ≈ zeros(T, B.dim, B.dim)
            end
        end
    end
end
