import numpy as np
from UtilsStaglineFastRun.StaglineFastRunSubsonic import *
import pdb
import sys

# ========================================================================
 
# ----------------------------------
# | INPUTS FOR STAGLINE SIMULATION |
# ----------------------------------
 
# Inputs
# ------
 
# Name of the mixture
mixture = "air_11"

# Number of species according to the mixture used
n_species = "11"
 
# Initial value of the CFL
cfl_val = "1d0"
 
# Temperature of the wall
Twall = "350"
 
# Method to change the CFL
cfl_adaptive = ".TRUE."
cfl_inter = ".FALSE."
Log_CFL = ".FALSE."

# Stop condition
residual = -2.0

# Restarting from previous simulation
restart = ".TRUE."

# Restart from previously converged air_7 simulation (Can be turned to True only if restart is True)
air_5_restart = ".FALSE."

# Management of CFL (air_7)
#CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  250,  500, 1000, 5000]
#Iter       = [   45,  95, 145, 195, 245, 295, 595, 795, 1245, 1445, 1745, 1845]

# Management of CFL (air_7) WITH RESTART FILES
#CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  250,  500, 1000, 5000]
#Iter       = [   45,  95, 145, 195, 245, 295, 395, 495,  595, 695,  795, 895]

# Management of CFL (air_11)
#CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  250,  500, 1000, 5000]
#Iter       = [   45,  95, 145, 195, 245, 395, 695, 995, 1545, 1745, 2045, 2145]

# Management of CFL (air_11)
CFL_range  = [0.001, 0.01,  0.1,  0.5,  1.0,  5.0, 10.0,    100,    500]
Iter       = [  395,  895, 1595, 2295, 4295, 6295, 8295,  10295,  14295]


# ========================================================================
 
# ----------------------------------------------
# | INITIAL CONDITIONS FOR STAGLINE SIMULATION |
# ----------------------------------------------

# Target Pressure [mbar]
target_pressure = [15, 200]

# Target Power [kW]
target_temp = [7000, 8000]

# Target dynamic pressure [Pa]
target_pdyn = [50, 400]

# ========================================================================

# Looping on the different conditions
for h, targ_p in enumerate(target_pressure):

    for j, temp in enumerate(target_temp):

        for k, pdyn in enumerate(target_pdyn): 
            
            # Temperature at the inlet [K]
            Tin = temp
            
            # Static pressure of the chamber [Pa]
            pc = targ_p * 100
            
            # Radius of the probe [m]
            R = 0.025
            
            # If velocity data is available, set these parameters accordingly.  
            # Otherwise, assign to "None".  
            uin = None  # [m/s]
            vin = None  # [m/s]
            
            # If dynamic pressure data is available, set this parameter accordingly.  
            # Otherwise, assign to "None".  
            pdyn = pdyn   # [Pa]

            # ========================================================================

            # -----------------------------------------------------------------------------
            # | Paths to the simulation, "Template_files" folders and Stagline executable |
            # -----------------------------------------------------------------------------
            
            # Simulation folder path
            stagline_simulations_path = f"/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/StaglineSim/{mixture}_sim/Pc={targ_p}"

            # Restart Simulation folder path
            stagline_restart_simulations_path = f"/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/StaglineSim/{mixture}_sim/Pc={targ_p}/*"
            
            # Input file template path
            input_template_path = "/home/jpe/VKI/Project/Stagline/Template_files"

            # Path to the stagline executable
            stagline_exe_path = "/home/jpe/VKI/Project/Stagline/bin/stagline"

            # Catalicity files path
            catalicity_files_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/StaglineSim/Catalicity_files"

            # ========================================================================

            # --------------------------------------------------
            # | Visualistaion plot for the residuals behaviour |
            # --------------------------------------------------

            # Plot visualisation
            res_plot_visu = True

            # ========================================================================

            # ----------------------------------
            # | Launch of Stagline simulatrion |
            # ----------------------------------
                            
            StaglineFastRunSubsonic(Tin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,stagline_restart_simulations_path,input_template_path,stagline_exe_path,catalicity_files_path,CFL_range,Iter,air_5_restart,res_plot_visu)








