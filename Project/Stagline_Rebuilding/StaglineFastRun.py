import numpy as np
from Utils_StaglineFastRun import *
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
cfl_val = "1.0d-1"
 
# Temperature of the wall
Twall = "350"
 
# Method to change the CFL
cfl_adaptive = ".TRUE."
cfl_inter = ".FALSE."
 
# ----------------------------------------------------------------------
 
# Initial conditions
# ------------------
 
# Temperature at the inlet
Txin = 10000
Tyin = 10000
 
# Static pressure of the chamber
pc = 5000
 
# Radius of the probe
R = 0.025
 
# If velocity data is available, set these parameters accordingly.  
# Otherwise, assign to "None".  
uin = None  
vin = None  
 
# If dynamic pressure data is available, set this parameter accordingly.  
# Otherwise, assign to "None".  
pdyn = 15
 
# ----------------------------------------------------------------------
 
# Paths to the "Simulations" and "Template_files" folders in "Stagline"
# ---------------------------------------------------------------------
 
# The user must create two separate folders inside the "Stagline" directory:
 
# 1) "Simulations" folder:
#    - Must be inside the "Stagline" directory.
#    - "Stagline" is the folder containing a "bin" folder with the Stagline executable and the "src" folder containing the source codes.
 
stagline_simulations_path = "/home/jpe/VKI/Project/Stagline/MixtureAnalysis"
 
# 2) "Template_files" folder:
#    - Must be inside the "Stagline" directory.
#    - Must contain example input files for "air_5" and "air_11".
#    - These files must be named:
#        - Example_input_air_5
#        - Example_input_air_11
#    - Must contain The mesh.dat file.
 
input_template_path = "/home/jpe/VKI/Project/Stagline/Template_files"
 
# ----------------------------------------------------------------------
 
# End of user configuration section
# /!\ No modifications are required beyond this point /!\
 
 
 
 
 
 
 
 
 
 
 
 
 
#############################################################################
#############################################################################
#############################################################################
#############################################################################
 
# -----------------------------------------------
# | START OF THE SCRIPT FOR STAGLINE SIMULATION |
# -----------------------------------------------
 
# Processing user inputs
# ----------------------
inputs,init_cond,sim_name = UserInputsProcessing(Txin,Tyin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive)
 
# Runing script for input file modification
# -----------------------------------------
sim_folder_path = InputFileGenerator(stagline_simulations_path,input_template_path,sim_name,inputs,init_cond)
 
# Stagline simulation
# -------------------
RunStagline(sim_folder_path)







