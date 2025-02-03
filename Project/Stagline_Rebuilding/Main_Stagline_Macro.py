import numpy as np
from Utils_Macro import *

# -------------------------------------------------------
# | INPUTS / INITIAL CONDITIONS FOR STAGLINE SIMULATION |
# -------------------------------------------------------

# Inputs
# ------
mixture = "air_5"
n_species = "5"
cfl_val = "1.0d-1"
Twall = "350"
cfl_adaptive = ".TRUE."
cfl_inter = ".FALSE."

#-----------------------------

# Initial conditions
# ------------------
uin = -225.68437042842623
vin = 225.68437042842623
Txin = 5779.0
Tyin = 5779.0

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
density = MPPDensities(Txin,101300,"air_5")

# Initial conditions
# ------------------
init_cond = {} # creating the inputs dictionnary

# Initial conditions values
init_cond["Densities"] = density
init_cond["Velocities"] = np.array([uin,vin])
init_cond["Temperatures"] = np.array([Txin,Tyin])


# Simulation name
# ---------------
sim_name = f"sim_T={Txin}_V={np.abs(uin):.2f}"
inputs["Simulation_Name"] = sim_name

# Runing script for input file modification
# -----------------------------------------
sim_folder_path = InputFileGenerator(sim_name,inputs,init_cond)

# Stagline simulation
# -------------------
RunStagline(sim_folder_path)






