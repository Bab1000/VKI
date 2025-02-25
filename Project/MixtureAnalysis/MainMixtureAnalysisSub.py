import numpy as np
from UtilsStaglineFastRun.StaglineFastRunSubsonic import *
import pdb
import sys
 
# -------------------------------------------------------
# | INPUTS / INITIAL CONDITIONS FOR STAGLINE SIMULATION |
# -------------------------------------------------------
 
# Inputs
# ------
 
# Name of the mixture
mixture = "air_7"

# Number of species according to the mixture used
n_species = "7"
 
# Initial value of the CFL
cfl_val = "1.0d-3"
 
# Temperature of the wall
Twall = "350"
 
# Method to change the CFL
cfl_adaptive = ".FALSE."
cfl_inter = ".TRUE."
Log_CFL = ".FALSE."

# Stop condition
residual = -2.0

# Restarting from previous simulation
restart = ".TRUE."

# Restart from previously converged air_5 simulation (True only if restart is True)
air_5_restart = ".TRUE."

# Management of CFL
CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100, 500 ]
Iter       = [   45,  95, 145, 195, 245, 295, 1095, 1495, 1995]

#CFL_range  = [0.01, 0.1,   1,  10, 100, 500, 1000]
#Iter       = [  40,  90, 140, 190, 240, 540, 740]

#CFL_range  = [ 0.01, 0.1, 1.0, 10.0,  50, 100, 1000]
#Iter       = [   40,  90, 140,  190, 390, 590, 790]


 
# ========================================================================
 
# Initial conditions
# ------------------
 
# Temperature at the inlet [K]
T_test = [8200]
 
# Static pressure of the chamber [Pa]
pc = 5000
 
# Radius of the probe [m]
R = 0.025
 
# If velocity data is available, set these parameters accordingly.  
# Otherwise, assign to "None".  
uin = None  # [m/s]
vin = None  # [m/s]
 
# If dynamic pressure data is available, set this parameter accordingly.  
# Otherwise, assign to "None".  
pdyn = 180   # [Pa]
 
# ========================================================================
 
# Paths to the "Simulations" and "Template_files" folders in "Stagline"
# ---------------------------------------------------------------------
 
# The user must create two separate folders inside the "Stagline" directory:
 
# 1) "Simulations" folder:
#    - Must be inside the "Stagline" directory.
#    - "Stagline" is the folder containing a "bin" folder with the Stagline executable and the "src" folder containing the source codes.
 
stagline_simulations_path = "/home/jpe/VKI/Project/Stagline/SimulationsMixtureAnalysis"
 
# 2) "Template_files" folder:
#    - Must be inside the "Stagline" directory.
#    - Must contain example input files for "air_5" and "air_11".
#    - These files must be named:
#        - Example_input_air_5
#        - Example_input_air_11
#    - Must contain The mesh.dat file.
 
input_template_path = "/home/jpe/VKI/Project/Stagline/Template_files"
 
# ========================================================================
 
# End of user configuration section
# /!\ No modifications are required beyond this point /!\
 
#############################################################################
#############################################################################
#############################################################################
#############################################################################

for T in T_test:

    StaglineFastRunSubsonic(T,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,input_template_path,CFL_range,Iter,air_5_restart)
    







