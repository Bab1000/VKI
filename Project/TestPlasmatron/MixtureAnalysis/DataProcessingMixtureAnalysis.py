import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Utils_DataProcessingMixtureAnalysis import *

# ---------------------------------------------------
# | DATA PROCESSING FOR STAGLINE MASS FLOW ANALYSIS |
# ---------------------------------------------------

# ========================================================================

# Path for data processing
# ------------------------

# Global path of the results
glob_path = "/home/jpe/VKI/Project/MixtureAnalysis/Simulations"

# ========================================================================

# Looping on every targeted pressure and power
# --------------------------------------------

# Target pressure for data analysis
target_pressure = 50

# Name pattern of the air_5 simulations data folder  
air_5_simulation_folders_pattern = f"*_mix=air_5_{target_pressure}mbar"

# Name pattern of the air_7 simulations data folder  
air_7_simulation_folders_pattern = f"*_mix=air_7_{target_pressure}mbar"

# General path for air_5 simulation folders
air_5_simulation_general_path = os.path.join(glob_path, air_5_simulation_folders_pattern)

# General path for air_7 simulation folders
air_7_simulation_general_path = os.path.join(glob_path, air_7_simulation_folders_pattern)

# Locate all air_5 simulation folders
list_path_air_5 = glob.glob(air_5_simulation_general_path, recursive=True)

# Locate all air_7 simulation folders
list_path_air_7 = glob.glob(air_7_simulation_general_path, recursive=True)

# Sorting the file with respect to the temperatures values  

list_path_air_5 = sorted(list_path_air_5,key=extract_T_value)

list_path_air_7 = sorted(list_path_air_7,key=extract_T_value)


# ========================================================================

# Gathering the u,v for each files
# --------------------------------

q_tot_air_5 = []
q_tot_air_7 = []

for air_5_folder_path in list_path_air_5:

    # Locate the surface.dat file
    surface_file_general_path = os.path.join(air_5_folder_path,"*_surface.dat")

    surface_file_path = glob.glob(surface_file_general_path, recursive=True)

    # extracting Heat Flux value from the .dat file
    q_tot = extract_q_tot(surface_file_path[0])

    # Storing the value
    q_tot_air_5.append(q_tot)

for air_7_folder_path in list_path_air_7:

    # Locate the surface.dat file
    surface_file_general_path = os.path.join(air_7_folder_path,"*_surface.dat")

    surface_file_path = glob.glob(surface_file_general_path, recursive=True)

    # extracting Heat Flux value from the .dat file
    q_tot = extract_q_tot(surface_file_path[0])

    # Storing the value
    q_tot_air_7.append(q_tot)

# Transforming U and V in arrays
q_tot_air_5 = np.array(q_tot_air_5)  
q_tot_air_7 = np.array(q_tot_air_7)  

# ========================================================================

# ------------------------------
# | DATA LOADING FROM CSV FILE |
# ------------------------------

# Path to the CSV
CSV_path = "/home/jpe/VKI/Project/MixtureAnalysis/tests_Justin.xlsx"

# Define the columns to check for NaN values
columns_to_check = ["Pressure[mbar]", "massflow [g/s]", "Power[kW]", "HeatFlux(HS50mm)[kW/m2]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)"]

# Data loading
pressure,massflow,power,heat_flux,pitot,temperature = CSVReader(CSV_path,columns_to_check)

# Range
range = 2

exp_HF_mdot10 = []
exp_HF_mdot16 = []
exp_HF_mdot20 = []

T_modt10 = []
T_modt16 = []
T_modt20 = []

# Check only for target pressure
for i,press in enumerate(pressure):
    if (press >= target_pressure - range and press <= target_pressure + range):
        if (massflow[i] >= 10 - range and massflow[i] <= 10 + range):
            exp_HF_mdot10.append(heat_flux[i])
            T_modt10.append(temperature[i])

        elif (massflow[i]>= 16 - range and massflow[i] <= 16 + range):
            exp_HF_mdot16.append(heat_flux[i])
            T_modt16.append(temperature[i])

        elif (massflow[i]>= 20 - range and massflow[i] <= 20 + range):
            exp_HF_mdot20.append(heat_flux[i])
            T_modt20.append(temperature[i])

exp_HF_mdot10 = np.array(exp_HF_mdot10)
exp_HF_mdot16 = np.array(exp_HF_mdot16)
exp_HF_mdot20 = np.array(exp_HF_mdot20)

T_modt10 = np.array(T_modt10)
T_modt16 = np.array(T_modt16)
T_modt20 = np.array(T_modt20)

# ========================================================================



# Plotting results
# ----------------

# Path for saving the figure
plot_path_HF = os.path.join(glob_path, "Results", f"HF_analysis.pdf")

# Path for saving the figure
plot_path_HF_diff = os.path.join(glob_path, "Results", f"HF_diff_analysis.pdf")

# Temperature at the inlet [K]
T = np.arange(5000, 10001, 500)

# Plotting beta
plt.figure()
plt.plot(T,np.abs(q_tot_air_5)/1000)
plt.plot(T,np.abs(q_tot_air_7)/1000)
plt.plot(T_modt10,exp_HF_mdot10,'*')
plt.plot(T_modt10,exp_HF_mdot16,'*')
plt.plot(T_modt10,exp_HF_mdot20,'*')
plt.legend(["HF air_5", "HF air_7"])
plt.xlabel("T [K]")
plt.ylabel("HF [kW/m2]")
plt.title(f"Mixture Analysis")
plt.grid(True)  
plt.savefig(plot_path_HF, format='pdf', dpi=600, bbox_inches='tight', transparent=True)

# Plotting beta
plt.figure()
plt.plot(T,np.abs(q_tot_air_5 - q_tot_air_7)/1000)
plt.xlabel("T [K]")
plt.ylabel("HF error [kW/m2]")
plt.title(f"Heat Flux error between air_7 and air_5")
plt.grid(True)  
plt.savefig(plot_path_HF_diff, format='pdf', dpi=600, bbox_inches='tight', transparent=True)












