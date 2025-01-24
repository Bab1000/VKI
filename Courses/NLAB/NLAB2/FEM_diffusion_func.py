import numpy as np
from Utils_FEM import *
import pdb

def FEM(n_elem,mesh,nodes_list,elem_list,source,K_glob,M_glob,F_glob):

    # Local form fucntions evaluation
    for i in range(n_elem):

        # Get the element
        elem = elem_list[i]
        
        # Get node 1
        node1 = nodes_list[elem[0]]

        # Get node 2
        node2 = nodes_list[elem[1]]

        # Get location of node 1
        x1 = mesh[node1]

        # Get location of node 2
        x2 = mesh[node2]

        # Element length
        e_length = x2-x1

        # Evaluate K_e 
        K_e = K_e_eval(e_length)

        # Evaluate K_e 
        M_e = M_e_eval(e_length)

        # Evaluate F_e
        source1 = source[node1]
        source2 = source[node2]
        F_e = F_e_eval(source1,source2,x1,x2,e_length)

        # Assembly of K_glob
        K_glob[i : i + len(K_e), i : i + len(K_e)] += K_e

        # Assembly of F_glob
        F_glob[i : i + len(F_e)] += F_e

        # Assembly of M_glob
        M_glob[i : i + len(M_e), i : i + len(M_e)] += M_e

    return K_glob,F_glob,M_glob

def K_e_eval(e_length):
        
        K11 = 1/e_length
        K12 = -1/e_length
        K21 = -1/e_length
        K22 = 1/e_length

        K_e = np.array([[K11,K12],[K21,K22]])

        return K_e

def F_e_eval(source1,source2,x1,x2,e_length):
    # Evaluate F_e
    F1 = -source1/e_length * ((x2**2/2 - x2**2) - (x1**2/2 - x1*x2))
    F2 = source2/e_length * ((x2**2/2 - x2*x1) - (x1**2/2 - x1**2))
    F_e = np.array([[F1],[F2]])

    return F_e

def M_e_eval(e_length):
    M11 = e_length / 3
    M12 = e_length / 6
    M21 = e_length / 6
    M22 = e_length / 3

    M_e = np.array([[M11,M12],[M21,M22]])

    return M_e

def Backward_euler(n_elem,LHS,Dirichlet,M_glob,nt,dt):

    sol = []

    #Initial condition
    T0 = np.full(n_elem+1, 300).reshape(-1,1)

    T0[0] = Dirichlet
    
    # compute the missing stiffness from applying Dirichlet BC
    Missing_stiff = LHS[1:,0]*Dirichlet

    # compute LHS of the equation
    LHS = LHS[1:,1:]

    for n in range(nt):

        RHS = Applying_DirichletBC_time_dependent(T0,M_glob,dt,Missing_stiff)

        new_T0 = np.linalg.solve(LHS,RHS)

        new_T0_full = np.concatenate([[Dirichlet],new_T0.flatten()])

        sol.append(new_T0_full)

        T0 = new_T0_full.reshape(-1,1)

    return sol

def Time_dependent_diffusion(n_elem,mesh,nodes_list,elem_list,nt,dt):
     
   # Create the global K matrix
    K_glob = np.zeros((len(nodes_list),len(nodes_list)))

    # Create the global K matrix
    M_glob = np.zeros((len(nodes_list),len(nodes_list)))
    
    # Create the global F vector
    F_glob = np.zeros((len(nodes_list),1))

    # Compute source term
    source = Source_term(mesh)

    print("   ---> Building matrices ...")

    K_glob,F_glob,M_glob = FEM(n_elem,mesh,nodes_list,elem_list,source,K_glob,M_glob,F_glob)

    K_glob = np.array(K_glob)
    F_glob = np.array(F_glob)
    M_glob = np.array(M_glob)

    print(f"   ---> Applying BCs :")
    # --> Dirichlet BC

    # compute LHS of the equation
    LHS = (M_glob/dt + K_glob)

    # Dirichlet conditions
    Dirichlet = 500

    print(f"        \033[34mDirichlet (left wall)\033[0m  : {Dirichlet} K")
    print(f"        \033[34mNeumann   (right wall)\033[0m : dt/dx = 0")

    # Backward Euler
    print("Simulating time dependent temperature diffusion")

    sol = Backward_euler(n_elem,LHS,Dirichlet,M_glob,nt,dt)

    print("End of simulation ...")

    return sol   

def Diffusion(n_elem,mesh,nodes_list,elem_list):

    # Create the global K matrix
    K_glob = np.zeros((len(nodes_list),len(nodes_list)))

    # Create the global K matrix
    M_glob = np.zeros((len(nodes_list),len(nodes_list)))
    
    # Create the global F vector
    F_glob = np.zeros((len(nodes_list),1))

    # Compute source term
    source = Source_term(mesh)

    print("   ---> Building matrices ...")

    K_glob,F_glob,M_glob = FEM(n_elem,mesh,nodes_list,elem_list,source,K_glob,M_glob,F_glob)

    # Applying BCs
    K_glob = np.array(K_glob)
    F_glob = np.array(F_glob)

    print(f"   ---> Applying BCs :")
    
    # --> Dirichlet
    T0 = 500
    
    Missing_stiff = K_glob[1:,0]*T0
    K_glob = K_glob[1:,1:]
    F_glob = F_glob[1:]
    F_glob = F_glob - Missing_stiff.reshape(-1,1)

    print(f"        \033[34mDirichlet (left wall)\033[0m  : {T0} K")
    print(f"        \033[34mNeumann   (right wall)\033[0m : dt/dx = 0")
    
    # -- >Neumann nothing to do !!  dT/dx * phi(L) = 0 

    # Resolve the equation
    print("Simulating temperature diffusion")

    sol = np.linalg.solve(K_glob,F_glob)

    sol = np.concatenate([[T0],sol.flatten()])

    print("End of simulation ...")

    return sol

def Applying_DirichletBC_time_dependent(T0,M_glob,dt,Missing_stiff):

    # compute RHS
    RHS = (M_glob / dt) @ T0
    RHS = RHS[1:]
    RHS = RHS - Missing_stiff.reshape(-1,1)

    return RHS









