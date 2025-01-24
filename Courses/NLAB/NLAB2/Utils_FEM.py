import numpy as np
import math as m
import matplotlib.pyplot as plt
import os

def OneD_uniform_mesh_generation(n_elem,l):

    # 1D uniform mesh
    n_nodes = n_elem + 1

    mesh = np.linspace(0,l,n_nodes)

    # list of nodes
    nodes_list = []

    for i in range(n_nodes):
        nodes_list.append(i)

    # list of elements (list with the nodes associated to this element)
    elem_list = []

    for i in range(n_elem):
        elem_list.append([nodes_list[i],nodes_list[i+1]])

    # list of DOF (here only one DOF per node ==> DOF_list = nodes_list)
    DOF_list = nodes_list.copy()

    return mesh,nodes_list,elem_list,DOF_list

def GaussQuadrature():
    xi=0
    w=2
    return xi,w

def GaussInverseMapping(x1,l):

    xi_InvMapped = l/2 * (x1 + 1)

    return xi_InvMapped

def Source_term(mesh):
    source = -4000/np.exp(10) * np.exp(10 * (mesh - 1)**2) * (1 + 20 * (mesh - 1)**2)

    return source

def plot_diffusion(mesh,res):
    print("Plotting temperature diffusion simulation results...")

    os.makedirs("Results/Diffusion", exist_ok=True)
    save_path = os.path.join("Results", "Diffusion", "Diffusion.pdf")

    T_exact = 200 * np.exp(10*(mesh - 1)**2)/np.exp(10) + 300

    plt.figure()
    plt.plot(mesh, T_exact, linestyle='-',color='red', linewidth=1.5, label='Exact Temperature Profile')
    plt.plot(mesh, res, linestyle='--', color='blue',linewidth=1.5, label='Temperature Profile (diffusion problem)')
    plt.xlabel('Mesh size', fontsize=12)
    plt.ylabel('Temperature [K]', fontsize=12)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in: {save_path}")

def plot_diffusion_time_dependent(t_vec,mesh,sol_TD):

    print("Plotting time dependent temperature diffusion simulation results...")

    os.makedirs("Results/TD_Diffusion", exist_ok=True)
    save_path = os.path.join("Results", "TD_Diffusion", "TD_Diffusion.pdf")

    plt.figure()
    for solution, t in zip(sol_TD, t_vec):

        plt.plot(mesh, solution, label=f't = {t + 0.1:.1f}s',linestyle='-',linewidth=1.5)

    plt.xlabel('Mesh size', fontsize=12)
    plt.ylabel('Temperature', fontsize=12)
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in: {save_path}")

def print_mesh(n_elements):

    if n_elements < 1:
        print("The number of elements must be at least 1.")
        return
    
    n_nodes = n_elements + 1
    
    if n_nodes <= 10:
        # Full display if there are 10 or fewer nodes
        mesh_line = "o---" * n_elements + "o"
        node_numbers = "0  " + "  ".join(f"{node+1:2}" for node in range(n_nodes-1))
    else:
        # Display first 3 and last 3 for larger meshes
        mesh_line = "o---o---o--- ... ---o---o---o"
        node_numbers = (
            "0   1   2          "
            + "  ".join(f"{node:2}" for node in range(n_nodes - 3, n_nodes))
        )
    
    # Print the mesh and node numbers
    print(mesh_line)
    print(node_numbers)



























