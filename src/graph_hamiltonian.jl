using BoseHubbardDiagonalize
using JeszenszkiBasis
using LightGraphs #might be ueful later for finding nearest neighbour on hex lattice
using GraphPlot
using SparseArrays: sparse
using DocStringExtensions

function graph_hamiltonian(basis::AbstractSzbasis, g::AbstractGraph)
    rows=Int64[]
    cols=Int64[]
    elements=Float64[]
    U=16
    t=4
    L = nv(g)
    
    

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
            
            nebs = neighbors(g, j)
            for j_next in nebs
                
           
              

                for (site1,site2) in [(j,j_next)]
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
    end

    sparse(rows,cols,elements,length(basis),length(basis))

end
export graph_hamiltonian







  







   



 

        
     



