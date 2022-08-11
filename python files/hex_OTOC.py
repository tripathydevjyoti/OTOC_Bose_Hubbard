from quspin.basis import boson_basis_1d #basis for Bose Hubbard
from quspin.operators import hamiltonian
import numpy as np
import networkx as nx #for generating a hexagonal grid
import matplotlib.pyplot as plt # plotting 
from quspin.tools.measurements import ED_state_vs_time #for time evolution of state vector
from math import*
import scipy
from numpy import linalg as la

def Hamil(N1):
    valid_state = []
   

    #To create larger lattice, increase the values of m and n
    #i.e to create a grid of 2 hexagons m=1;n=2

    m =  1 # number of rows of hexagons in the lattice
    n =  1 # number of columns of hexagons in the lattice
    isPBC = False # if True, use periodic boundary conditions
#
### build graph using networkx
    hex_graph = nx.generators.lattice.hexagonal_lattice_graph(m, n, periodic=isPBC)
# label graph nodes by consecutive integers
    hex_graph = nx.convert_node_labels_to_integers(hex_graph)
# set number of lattice sites
    N = hex_graph.number_of_nodes()
    print('constructed hexagonal lattice with {0:d} sites.\n'.format(N))
# visualise graph
    pos = nx.spring_layout(hex_graph, seed=42, iterations=100)
    nx.draw(hex_graph, pos=pos, with_labels=True)
    plt.show(block=False)
    
    L=N #number of sites in the lattice is equal to number of nodes in the graph
    
    #U/J=4
    t = 4 # hopping coefficient
    U = 16 # on-site fermion interaction strength
    mu = 0 #no chemical potential in the hamiltonian
##### set up Bose-Hubbard Hubbard Hamiltonian with quspin #####

    # creates a boson basis for L sites with N1 no of bosons.
    basis = boson_basis_1d(L,Nb=N1,sps=N1+1)  #states per site(sps) also includes having 0 boson state. hence the N1+1
    print('Hilbert space size: {0:d}.\n'.format(basis.Ns))


    # site-coupling lists
    hop   = [[-t, i, j] for i in range(N) for j in hex_graph.adj[i]] #for a 1d lattice j would be i+1
    #hop   = [[-t,0,1], [-t,1,2], [-t,2,3], [-t,3,4], [-t,4,5] ]
    #hop=[[-t,i,(i+1)%L] for i in range(L)] #PBC
    interact = [[0.5* U, i, i] for i in range(N)]
    pot=[[-mu-0.5*U,i] for i in range(L)]
    
    # define site-coupling lists [hermitian conjugates "-+|" and "|-+" contained in tunnelling list]
    
    #the ['-+',hop] term is in comments to avoid double counting when j runs over n.n of each vertex.
    static=[['+-',hop],['nn',interact],['n',pot]]#,['-+',hop]]
    dynamic=[] #no time dependent term in the hamiltonian
#
    #constructing Hamiltonian  
    temp_mat=hamiltonian(static,dynamic,basis=basis,dtype=np.float64)
    H=temp_mat.todense()
    
    #takes the basis generated , stores and returns them as lists of integers for easier creation and annhilation
    for i in range(basis.Ns):
        s = basis[i]
        s_str = basis.int_to_state(s)
        b_str = ""
        for i in range(len(s_str)):
            if i!=0 and i!=(len(s_str)-1):
                b_str=b_str+s_str[i]
        b = []
    
        for i in range(len(b_str)):
            if b_str[i]!=" ":
                b.append(int(b_str[i]))
        
        valid_state.append(b) 

    return valid_state,H    



def create(ket,i):
    ls = [ket[l] for l in range(len(ket))]
    ls[i] = ket[i]+1
   
    return(ls)
    
    

def destroy(ket,i):
    
    #Note that if site i has 0 bosons, this action will give you -1 after operation.
    #Take care to put at least 1 boson on the annihilation site 
    ls = [ket[l] for l in range(len(ket))]
    ls[i] = ket[i]-1
    
    return(ls)


#Writing this code for one hexagon with 6 bosons.

nb=6 #no of bosons


list_states,H_mat=Hamil(nb)
list_states1,H_mat1=Hamil(nb-1)#used to evolve states after annihilating once
list_states2,H_mat2=Hamil(nb-2)#used to evolve states after annihilating twice
   
    

#new and new2 used for some storage   
new=np.zeros(len(list_states2),dtype='complex') 
new2=np.zeros(len(list_states1),dtype='complex')
new3=np.zeros(len(list_states),dtype='complex')


#computes all the eigenvalues required for time evolution in the otoc function
E,V=la.eigh(H_mat)
E1,V1=la.eigh(H_mat1)
E2,V2=la.eigh(H_mat2)

def OTOC(dsj,time):


    init_state=[1,1,1,1,1,1] #change this inital state according to the size of lattice and no of bosons.
    
    #these three matrices are created to draw unit vectors and make superpositions after time evolution.
    track=np.identity(len(list_states),dtype='complex')
    track1=np.identity(len(list_states1),dtype='complex')
    track2=np.identity(len(list_states2),dtype='complex')
    
    #find position of initial state among basis vectors
    for i in range(len(list_states)):
        if(str(init_state)==str(list_states[i])):
            a=i
    
    #initialise the state vector based on the index found
    init_vec1=track[:,a]
    
    #apply annihilation op at i'th site
    new_state=destroy(init_state,dsj) 


    # see which of the new basis vectors does the annhilated state correspond to
    for i in range(len(list_states1)):
        if(str(new_state)==str(list_states1[i])): 
            a=i
            
        
    #initialise state vector for first time evolution
    init_vec=track1[a,:]
  
    
    
    t_list=[time] #forward evolution
    
    #ED_state_vs_time is a quspin fucntion that evolves a given initial state with the help of eigenvalues and eigenvectors
    new_vec=ED_state_vs_time(init_vec,E1,V1,t_list,iterate=False)
    new_vec=new_vec.tolist()
   
  
  
    
    
    store=[] #list for storing the superposition information
    
    
   
    
    for i in range(len(list_states1)):
        
        if(np.vdot(track1[:,i],new_vec)>0.001):
            store.append([i,np.vdot(track1[:,i],new_vec)]) #stores index of full basis vec and its coeff in superposition expression
    
    number=0
    
    temp_vec=np.zeros(len(track1[:,a])) #used to store inital state for next time evolution. but before that, there is another annihilation operator
    
    for j in range(len(store)):
        # the previous evolution results in a superposition of basis vectors. we then act annhilation operator on each of these
        # states , find their index and create the final state by adding them with appropriate coeff appearing in the superposition.
        temp_state=destroy(list_states1[store[j][0]],0)
      
        
        if(np.sum(temp_state)>0):
            number=number+1
#            
            
            for i in range(len(list_states2)):
                if(str(temp_state)==str(list_states2[i])):
                    a=i
        
        
        if(number==1):
            for i in range(len(track2[:,a])):
                track2[:,a][i]=store[j][1]*track2[:,a][i]
                
            temp_vec=track2[:,a]
            
            
        elif(number>1):
            
            for i in range(len(track2[:,a])):
                track2[:,a][i]=store[j][1]*track2[:,a][i]
            temp_vec=temp_vec+track2[:,a]   
            
            
        
        
    t_list=[-time]  #backward evolution 
    new_vec1=ED_state_vs_time(temp_vec,E2,V2,t_list,iterate=False) #another time evolution
    
    
    
    #from here on, we will be doing the same thing we did for annihilation operators but this time with creation operators.
    store1=[]
    
    for i in range(len(new_vec1)):
        new[i]=new_vec1[i][0]

                   
    
    for i in range(len(list_states2)):
        
        if(np.vdot(track2[:,i],new))>0.001:
            store1.append([i,np.vdot(track2[:,i],new)])
    
   
    number=0  
    temp_vec1=np.zeros(len(track1[:,a]))
    for j in range(len(store1)):
        
        temp_state=create(list_states2[store1[j][0]],dsj)
        
        
        if(np.sum(temp_state)>0):
            number=number+1
#             
            
            for i in range(len(list_states1)):
                if(str(temp_state)==str(list_states1[i])):
                    
                    a=i
                 
        if(number==1):
            for i in range(len(track1[:,a])):
                track1[:,a][i]=store1[j][1]*track1[:,a][i]
                
            temp_vec1=track1[:,a]
            
            
        elif(number>1):
            
            for i in range(len(track1[:,a])):
                track1[:,a][i]=store1[j][1]*track1[:,a][i]
            temp_vec1=temp_vec1+track1[:,a]   
            
        

    
    
    t_list=[time]
    new_vec2=ED_state_vs_time(temp_vec1,E1,V1,t_list,iterate=False)
    
    
    store2=[]
    for i in range(len(new_vec2)):
        new2[i]=new_vec2[i][0]
    
    for i in range(len(list_states1)):
        if(np.vdot(track1[:,i],new2))>0.001:
            store2.append([i,np.vdot(track1[:,i],new2)])
    temp_vec2=np.zeros(len(track[:,a]))         
    for j in range(len(store2)):
        temp_state=create(list_states1[store2[j][0]],0)
        if(np.sum(temp_state)>0):
            number=number+1
            for i in range(len(list_states)):
                if(str(temp_state)==str(list_states[i])):
                    
                    a=i
              
        if(number==1):
            for i in range(len(track[:,a])):
                track[:,a][i]=store2[j][1]*track[:,a][i]
                
            temp_vec2=track[:,a]
            
            
        elif(number>1):
            
            for i in range(len(track[:,a])):
                track[:,a][i]=store2[j][1]*track[:,a][i]
            temp_vec2=temp_vec2+track[:,a]   
        
    t_list=[-time]  
    new_vec3=ED_state_vs_time(temp_vec2,E,V,t_list,iterate=False)
    for i in range(len(new_vec3)):
        new3[i]=new_vec3[i][0]
    return ((abs(np.vdot(init_vec1,new3))) )
    

OTOC_list4=[]
t_list1=np.linspace(0,0.2,100)
#t_list1=np.linspace(0,0.1,2)
t_list_new=[4*x for x in t_list1]
for t in t_list1:
    
    
    #OTOC_list2.append(OTOC(2,t))
    #OTOC_list3.append(OTOC(3,t))
    #OTOC_list4.append(OTOC(4,t))
    OTOC_list4.append(OTOC(4,t))
    #print(t)
    #OTOC_list6.append(OTOC(6,t))
    

    
plt.plot(t_list_new,OTOC_list4)
plt.xlabel('Jt')
plt.ylabel('OTOC')
#plt.scatter(t_list_new,OTOC_list5)
plt.grid()

    