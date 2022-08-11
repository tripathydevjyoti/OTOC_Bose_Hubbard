using BoseHubbardDiagonalize
using JeszenszkiBasis
using LightGraphs #might be ueful later for finding nearest neighbour on hex lattice
using GraphPlot
using SparseArrays: sparse

function sparse_hamiltonian(basis::AbstractSzbasis)
    rows=Int64[]
    cols=Int64[]
    elements=Float64[]

    for (i,bra) in enumerate(basis)
        #diagonal entries
        Usum=0
        for j in 1:basis.K
            Usum += bra[j]*(bra[j]-1)
        end
        push!(rows,i)
        push!(cols,i)
        push!(elements,U*Usum/2)

        #off-diagonal entries

        for j in 1:L
            j_next = mod(j+1, 1:basis.K)

            for (site1,site2) in [(j,j_next),(j_next,j)]
                if bra[site1]>0
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis
                        push!(rows, i)
                        push!(cols, serial_num(basis,ket))
                        push!(elements, -t*sqrt(bra[site1])*sqrt(bra[site2]+1))
                    end
                end
            end
        end
    end

    sparse(rows,cols,elements,length(basis),length(basis))

end


L = 6 #number of sites in the lattice
N = 6 #number of bosons in the lattice
basis1=Szbasis(L,N)
basis2=Szbasis(L,N-1)
basis3=Szbasis(L,N-2)

const U=16
const t=4
H1 = sparse_hamiltonian(basis1)
H2 = sparse_hamiltonian(basis2)
H3 = sparse_hamiltonian(basis3)



function time_evoultion(eig_vals,eig_vecs,init_state,time)
    coeff_list = Float64[]
    out_vec = zeros(length(init_state))
    for i in 1:length(eig_vals)
        conj_vec = conj!(eig_vec[i])
        coeff = exp(-1im*eig_vals[i]*time)* ( vecdot(conj_vec,init_state) )
        out_vec = out_vec + coeff*eig_vecs[i]
    end
    return out_vec

end

function superposition(inp_vec)
    ket_list=[]
    coeff_list=[]
    track = zeros(len(inp_vec))
    for i in 1:len(inp_vec)
        track[i]=1
        if vecdot(track,inp_vec)>0.001
            push!(ket_list,basis[i])
            push!(coeff_list,vecdot(inp_vec))
        end
        track[i]=0
    end
    return ket_list,coeff_list
end             



function create(ket,i)
    ket[i] = ket[i]+1
    return ket
end

function destroy(ket,i)
    ket[i] = ket[i]-1
    return ket
end

function find_index(inp_state,basis)
    for i in 1:length(basis)
        if inp_state == basis[i]
            index =i
        end
    end
    return index 
end           


function OTOC(site1,site2,time)
    init_state1=[1,1,1,1,1,1]
    
    init_state = destroy(init_state1,site1)

    for i in 1:len(basis2)
        if init_state == basis2[i]
            index = i
        end
    end
    
    track = zeros(length(basis2))
    track[index]=1
    
    new_state = time_evoultion(E2,V2,track,time)
    
    superpos_coeff,superpos_state = superposition(new_state)

    new_state = zeros(length(basis3))
    for i in 1:length(superpos_state)
        
        superpos_state[i] = destroy[superpos_state[i],site2]

        track = zeros(length(basis3))
        index = find_index(superpos_state[i],basis3)
        track[index]=1
        new_state = new_state + superpos_coeff[i]*track
    end

    new_state = time_evoultion(E3,V3,new_state,-time)
    
    superpos_coeff,superpos_state = superposition(new_state)

    new_state = zeros(length(basis2))
    for i in 1:length(superpos_state)
        
        superpos_state[i] = create[superpos_state[i],site1]

        track = zeros(length(basis2))
        index = find_index(superpos_state[i],basis2)
        track[index]=1
        new_state = new_state + superpos_coeff[i]*track
    end

    new_state = time_evolution(E2,V2,new_state,time)

    superpos_coeff,superpos_state = superposition(new_state)

    new_state = zeros(length(basis1))
    for i in 1:length(superpos_state)
        
        superpos_state[i] = create[superpos_state[i],site2]

        track = zeros(length(basis1))
        index = find_index(superpos_state[i],basis1)
        track[index]=1
        new_state = new_state + superpos_coeff[i]*track
    end

    fin_state = time_evolution(E1,V1,new_state,-time)

    track=zeros(length(basis1))
    index = find_index(init_state1,basis1)
    track[index] = 1
    return vecdot(track,fin_state)

end







    
    






    

