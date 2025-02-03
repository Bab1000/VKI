import numpy as np
from Utils_Macro import *

# -------------------------------------------------------
# | INPUTS / INITIAL CONDITIONS FOR STAGLINE SIMULATION |
# -------------------------------------------------------

# Inputs
# ------
mixture = "air_11"
n_species = "11"
cfl_val = "1.0d-2"
Twall = "350"
cfl_adaptive = ".TRUE."
cfl_inter = ".FALSE."

#-----------------------------

# Initial conditions
# ------------------
uin = -377.15
vin = 377.15
Txin = 7237
Tyin = 7237
pc = 5000

#-----------------------------

# En of user input section of the scripts

#############################################################################
#############################################################################
#############################################################################
#############################################################################

# Inputs
# ------
inputs = {}  # creating the inputs dictionnary

# Inputs values
inputs["Mixture"] = mixture
inputs["Number_of_species"] = n_species
inputs["CFL_number"] = cfl_val
inputs["Twall"] = Twall
inputs["Inter_CFL"] = cfl_inter
inputs["Adaptive_CFL"] = cfl_adaptive

# Gathering the species density
# -----------------------------
density = MPPDensities(Txin,pc,mixture)

# Initial conditions
# ------------------
init_cond = {} # creating the inputs dictionnary

# Initial conditions values
init_cond["Densities"] = density
init_cond["Velocities"] = np.array([uin,vin])
init_cond["Temperatures"] = np.array([Txin,Tyin])


# Simulation name
# ---------------
sim_name = f"sim_T={Txin}_V={np.abs(uin):.2f}_{mixture}"
inputs["Simulation_Name"] = sim_name

# Runing script for input file modification
# -----------------------------------------
sim_folder_path = InputFileGenerator(sim_name,inputs,init_cond)

# Stagline simulation
# -------------------
RunStagline(sim_folder_path)






