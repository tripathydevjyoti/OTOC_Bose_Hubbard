M = N = 3
D = Int(factorial(N + M − 1) / factorial(N) / factorial(M − 1))

B = Basis(M, N)

@testset "Dense Hamiltonian" begin
    T = Float64
    J = T(1)
    U = T(1/2)
    H = hamiltonian(M, N, J, U, :OBC)

    @test size(H) == (D, D)
    @test transpose(H) == H
end
