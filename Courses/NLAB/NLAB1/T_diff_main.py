import numpy as np
from T_diff_func import *

# -----------------
# | GLOBAL INPUTS |
# -----------------

L = 1                    # length of the mesh
n_points = 160           # number ofpoints in the mesh
k = 225                  # diffusion constant
stencil = [-3,-2,-1,0,1] # spatial discretization stencil

#############################################################################
#############################################################################

# -------------------------------
# | SPATIAL SCHEME COEFFICIENTS |
# -------------------------------

print("\033[32mComputing the coefficients for the spatial schemes and k\033[0m")

# coefficient for the spatial derivative distretizations
coeff_Txx = FDM_coeff(stencil,derivative_order = 2)
coeff_Tx = FDM_coeff(stencil,derivative_order = 1)
coeff_k = FDM_coeff(stencil,derivative_order = 1)

print("   ---> Txx coefficients : " + str(coeff_Txx))
print("   ---> Tx coefficients  : " + str(coeff_Tx))
print("")

#############################################################################
#############################################################################

# -----------------------------
# | MESH CONVERGENCE ANALYSIS |
# -----------------------------

print("\033[32mMesh convergence analysis with periodic BC\033[0m")
print("Checking for solution convergence ...")

# Inputs for mesh convergence study 
beta = 0.2             # Neumann number
nt = 70001           # simulation time

# Periodic BC analysis
test_mesh = [10,20,40,80,160]
err,res_vec = Mesh_convergence(test_mesh,stencil,coeff_Txx,k,beta,nt)
convergence_slope = np.polyfit(np.log(test_mesh), np.log(err), 1)[0]

print("Computing convergence slope ...")
print("   ---> Convergence slope : " + str(convergence_slope))
print("")


#############################################################################
#############################################################################

# ----------------------
# | STABILITY ANALYSIS |
# ----------------------

print("\033[32mStability analysis with periodic BC\033[0m")
print("Checking for stability by varying beta ...")

# Inputs for mesh convergence study 
test_beta = np.arange(0.1, 0.51, 0.02)      # Neumann number
       
nt = 70001                                  # simulation time

# Periodic BC analysis
mesh = Mesh_generation(L,n_points)

residual_stability = Stability_analysis(test_beta,mesh,stencil,coeff_Txx,k,nt)
print("")

#############################################################################
#############################################################################

# -----------------------------
# | NON-PERIODIC BC ANALYSIS |
# -----------------------------

print("\033[32mTemperature diffusion analysis with Neumann and Dirichlet BC\033[0m")

# Inputs for non-periodic BC analysis
beta = 0.2                # Neumann number
nt = 70000                # simulation time

# mesh generation
mesh = Mesh_generation(L,n_points)

# computing initial condition
u0 = np.full(len(mesh),300)

# diffusion analysis
dx = mesh[1] - mesh[0]
dt =2*beta*dx**2/k
res = Temperature_diffusion_nonperiodicBC(u0,stencil,mesh,coeff_Txx,coeff_Tx,k,dx,dt,nt,dirichlet_bc=800)
print("")
print("End of simulation !")
print("")

#############################################################################
#############################################################################

# ---------------------------------------------------------------
# | MESH CONVERGENCE ANALYSIS ANALYSIS WITH K VARIABLE IN SPACE |
# ---------------------------------------------------------------

print("\033[32mMesh convergence analysis with periodic BC and variable k\033[0m")

# Inputs for mesh convergence study 
beta = 0.2             # Neumann number
nt = 70001           # simulation time

# Periodic BC analysis
test_mesh_k = [10,20,40,80,160]
err_k,residual_vec = Mesh_convergence_variable_k(test_mesh_k,stencil,coeff_Txx,coeff_Tx,coeff_k,k,beta,nt)
convergence_slope_k = np.polyfit(np.log(test_mesh), np.log(err_k), 1)[0]

print(np.shape(residual_vec))

print("Computing convergence slope ...")
print("   ---> Convergence slope : " + str(convergence_slope_k))
print("")
print("End of simulation ! ")
print("")

#############################################################################
#############################################################################

# ---------
# | PLOTS |
# ---------

print("\033[32mResults analysis\033[0m")


plot_mesh_convergence(err,test_mesh,convergence_slope,res_vec[-1])
plot_stability(residual_stability)
plot_nonperiodicBC(res,mesh,dt)
plot_mesh_convergence_variable_k(err_k,test_mesh_k,convergence_slope_k,residual_vec)
print("")
print("End of the script !")