using .OTOC_example
using Test
L =2
N =3
t = 1.0

ψ₀ = [1,2]

@testset "Compare OTOC function with analytics'" begin
    

    basis1 = Szbasis(L,N)
    basis2 = Szbasis(L,N-1)
    basis3 = Szbasis(L,N-2)


    H1 = sparse_hamiltonian(basis1, L)
    H2 = sparse_hamiltonian(basis2, L)
    H3 = sparse_hamiltonian(basis3, L)
    
    destroy(ψ₀,2)
    
    m2 = exp(-1im*Matrix(H2))
    v =(m2*(sqrt(2)*[0,1,0]))

    m3 = exp(+1im*Matrix(H3))
    v = m3*[sqrt(2)*(-0.7254232071350418 + 0.10668441034775135im),
        -0.6307931351627424 - 0.7258506973950558im]

    v = [ 0.0 + 0.0im,
        -0.5688574029231347 + 0.6021281395147168im,
        0.34089532057651484 + 1.584763798771675im]


    m4 = exp(-1im*Matrix(H2))

    v = m4*v

    v= [sqrt(3)*(1.1765821188141055 + 0.47369656372770386im),
       sqrt(2)*(0.2683636450488768 - 0.7637948073712356im),
       0.39386172946903375 - 0.9458226188097126im,0]



    m5 = exp(+1im*Matrix(H1))
    
    @test isapprox(abs((m5*v)[3]), OTOC_lattice([1,2],2,1,1.0,L,N))

end

