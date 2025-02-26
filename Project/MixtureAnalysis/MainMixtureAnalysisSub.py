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
mixture = "air_5"

# Number of species according to the mixture used
n_species = "5"
 
# Initial value of the CFL
cfl_val = "1.0d-3"
 
# Temperature of the wall
Twall = "350"
 
# Method to change the CFL
cfl_adaptive = ".TRUE."
cfl_inter = ".FALSE."
Log_CFL = ".FALSE."

# Stop condition
residual = -2.0

# Restarting from previous simulation
restart = ".FALSE."

# Restart from previously converged air_5 simulation (True only if restart is True)
air_5_restart = ".FALSE."

# Management of CFL
CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  200,  500,  750, 1000, 5000, 10000]
Iter       = [   45,  95, 145, 195, 245, 295, 595, 995, 2000, 2345, 2545, 2745, 2845,  2945]

#CFL_range  = [0.01, 0.1,   1,  10, 100, 500, 1000]
#Iter       = [  40,  90, 140, 190, 240, 540, 740]

#CFL_range  = [ 0.01, 0.1, 1.0, 10.0,  50, 100, 1000]
#Iter       = [   40,  90, 140,  190, 390, 590, 790]


 
# ========================================================================
 
# Initial conditions
# ------------------
 
# Temperature at the inlet [K]
T_test = np.arange(5000, 10001, 500)
 
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
pdyn = 100   # [Pa]
 
# ========================================================================

# -----------------------------------------------------------------------------
# | Paths to the simulation, "Template_files" folders and Stagline executable |
# -----------------------------------------------------------------------------

# Simulation folder path
stagline_simulations_path = "/home/jpe/VKI/Project/MixtureAnalysis/Simulations"

# Input file template path
input_template_path = "/home/jpe/VKI/Project/Stagline/Template_files"

# Catalicity files path
catalicity_files_path = "/home/jpe/VKI/Project/Stagline/Catalicity_files"

# Path to the stagline executable
stagline_exe_path = "/home/jpe/VKI/Project/Stagline/bin/stagline"

# ========================================================================

# --------------------------------------------------
# | Visualistaion plot for the residuals behaviour |
# --------------------------------------------------

# Plot visualisation
res_plot_visu = True

 
# End of user configuration section
# /!\ No modifications are required beyond this point /!\
 
#############################################################################
#############################################################################
#############################################################################
#############################################################################

for T in T_test:

    StaglineFastRunSubsonic(T,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,input_template_path,stagline_exe_path,catalicity_files_path,CFL_range,Iter,air_5_restart,res_plot_visu)
    







