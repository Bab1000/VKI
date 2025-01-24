import numpy as np
from FEM_diffusion_func import *
import matplotlib.pyplot as plt

# -----------------
# | GLOBAL INPUTS |
# -----------------

print("\033[32mInitialising simulation\033[0m")

n_elem = 100         # number of elements
l = 1                 # length of the mesh
order = 1             # order of the polynomials used as basis functions
t_end = 1            # Final time
dt = 0.1                # Time step size
nt = int(t_end / dt)  # Number of time steps

print(f"   ---> \033[34mNumber of elements in the mesh\033[0m : {n_elem}")
print(f"   ---> \033[34mLength of the domain\033[0m : {l}")
print(f"   ---> \033[34mTotal time of simulation\033[0m : {t_end}")
print(f"   ---> \033[34mTimestep size\033[0m : {dt}")
print(f"   ---> \033[34mTotal number of time steps\033[0m : {dt}")
print("Mesh type :")
print_mesh(n_elem)
print("")

#############################################################################
#############################################################################

# ------------------------------------
# | TEMPERATURE DIFFUSION SIMULATION |
# ------------------------------------

print("\033[32mComputing the temperature diffusion\033[0m")

# Mesh generation
print(f"Generating the mesh with {n_elem} elements ...")
mesh,nodes_list,elem_list,DOF_list = OneD_uniform_mesh_generation(n_elem,l)
print(f"   ---> \033[34mMesh {n_elem}\033[0m  successfully generated")

print("Simulating Temperature diffusion ...")
sol = Diffusion(n_elem,mesh,nodes_list,elem_list)

plot_diffusion(mesh,sol)
print("")

#############################################################################
#############################################################################

# ---------------------------------------------------
# | TIME DEPENDENT TEMPERATURE DIFFUSION SIMULATION |
# ---------------------------------------------------

print("\033[32mComputing the time dependent temperature diffusion\033[0m")

print("Simulating Temperature diffusion ...")
sol_TD = Time_dependent_diffusion(n_elem,mesh,nodes_list,elem_list,nt,dt)

t_vec = np.zeros(nt)

for i in range(nt):
    t_vec[i] = dt * i

plot_diffusion_time_dependent(t_vec,mesh,sol_TD)




















