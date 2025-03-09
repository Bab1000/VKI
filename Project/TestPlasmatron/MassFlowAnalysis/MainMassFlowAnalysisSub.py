import numpy as np
from UtilsStaglineFastRun.StaglineFastRunSubsonic import *
from Utils_MassFlowAnalysis import *
import pdb
import sys

# ========================================================================
 
# ----------------------------------
# | INPUTS FOR STAGLINE SIMULATION |
# ----------------------------------
 
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

# Restart from previously converged air_7 simulation (Can be turned to True only if restart is True)
air_5_restart = ".FALSE."

# Management of CFL (air_7 15mbar)
#CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  250,  500, 1000, 5000]
#Iter       = [   45,  95, 145, 195, 245, 295, 595, 795, 1245, 1445, 1745, 1845]

# Management of CFL (air_11 15mbar)
CFL_range  = [ 0.01, 0.1, 0.5, 1.0, 5.0,  10,  50, 100,  250,  500, 1000, 5000]
Iter       = [   45,  95, 145, 195, 245, 595, 795, 1245, 1445, 1745, 1845, 1995]
 
# ========================================================================

# ------------------------------
# | DATA LOADING FROM CSV FILE |
# ------------------------------

# Path to the CSV
CSV_path = "/home/jpe/VKI/Project/MassFlowAnalysis/tests_Justin.xlsx"

# Define the columns to check for NaN values
columns_to_check = ["Pressure[mbar]", "massflow [g/s]", "Power[kW]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)"]

# Data loading
pressure,massflow,power,pitot,temperature = CSVReader(CSV_path,columns_to_check)

# ========================================================================
 
# ----------------------------------------------
# | INITIAL CONDITIONS FOR STAGLINE SIMULATION |
# ----------------------------------------------

# Target Pressure [mbar]
target_pressure = [200]

# Target Power [kW]
target_power = [300,350]

# Tolerance
tolerance = 3.6

# ========================================================================

# Looping on the different conditions
for h, targ_p in enumerate(target_pressure):

    for i, targ_P in enumerate(target_power):  # Unpacking enumerate properly

        for j, pres in enumerate(pressure):  # Unpacking enumerate properly

            if (pres >= targ_p - tolerance and pres <= targ_p + tolerance and
                power[j] >= targ_P - tolerance and power[j] <= targ_P + tolerance):
    
                # Temperature at the inlet [K]
                Tin = temperature[j]
                
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
                pdyn = pitot[j]   # [Pa]

                # Mass flow [g/s]
                mdot = massflow[j]

                # Power [kW]
                Pw = power[j]
    
                # ========================================================================

                # -----------------------------------------------------------------------------
                # | Paths to the simulation, "Template_files" folders and Stagline executable |
                # -----------------------------------------------------------------------------
                
                # Simulation folder path
                stagline_simulations_path = f"/home/jpe/VKI/Project/MassFlowAnalysis/{mixture}_sim/Pc={targ_p}_Pw={targ_P}/mdot={mdot}"
                
                # Input file template path
                input_template_path = "/home/jpe/VKI/Project/Stagline/Template_files"

                # Path to the stagline executable
                stagline_exe_path = "/home/jpe/VKI/Project/Stagline/bin/stagline"

                # Catalicity files path
                catalicity_files_path = "/home/jpe/VKI/Project/Stagline/Catalicity_files"

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
                                
                StaglineFastRunSubsonic(Tin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,input_template_path,stagline_exe_path,catalicity_files_path,CFL_range,Iter,air_5_restart,res_plot_visu)
    







