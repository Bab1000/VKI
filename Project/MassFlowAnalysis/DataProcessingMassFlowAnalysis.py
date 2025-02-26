import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# ---------------------------------------------------
# | DATA PROCESSING FOR STAGLINE MASS FLOW ANALYSIS |
# ---------------------------------------------------

# ========================================================================

# Target pressure and power for the mass flow analysis
# ----------------------------------------------------

# Target pressure
targ_p = [50,100]

# Target power 
targ_P = [150,200,250,300,350]

# ========================================================================

# Path for data processing
# ------------------------

# Global path of the results
glob_path = "/home/jpe/VKI/Project/MassFlowAnalysis/air_5_sim"

# ========================================================================

# Looping on every targeted pressure and power
# --------------------------------------------
for press in targ_p:
    for i, Pow in enumerate(targ_P):

        # Name pattern of the data folder  
        pres_pow_pattern_name = f"Pc={press}_Pw={Pow}"

        # Name pattern of the mass flow folders
        massflow_pattern_name = "mdot=*"

        # Name pattern of the flowfield data file
        flowfield_pattern_name = "*flowfield.dat"

        # Full path to the pressure-power folder
        pres_pow_folder_path = os.path.join(glob_path, pres_pow_pattern_name)

        # Updated general path to match deeper folder structure
        flowfiled_general_path = os.path.join(pres_pow_folder_path, massflow_pattern_name, "**", flowfield_pattern_name)

        # Locate the flowfield data files
        list_path = glob.glob(flowfiled_general_path, recursive=True)

        # Sorting the file with respect to the massflow values
        def extract_mdot_value(filepath):
            match = re.search(r"mdot=([\d\.]+)", filepath)  
            return float(match.group(1)) if match else float('inf')  

        list_path = sorted(list_path,key=extract_mdot_value)
        


        # ========================================================================

        # Gathering the u,v for each files
        # --------------------------------

        # Initialisation of the storing matrices
        U = []
        V = []

        for i, flowfield_path in enumerate(list_path):

            # Reading the file
            df = pd.read_csv(flowfield_path, delim_whitespace=True, header=None,skiprows=4)

            # Extracting u and v
            u = df.iloc[:, 1].tolist()
            v = df.iloc[:, 2].tolist()
            mesh = df.iloc[:, 0].tolist()

            # Store the values in a matrix
            U.append(u)
            V.append(v)

        # Transforming U and V in arrays
        U = np.array(U)  
        V = np.array(V)  
        mesh = np.array(mesh)

        # ========================================================================

        # Computation of the velocity gradient
        # ------------------------------------
        # Probe radius [m]
        R = 0.025

        # Computing the velocity gradient
        beta_10_g_s = np.abs((U[0][:] + V[0][:])/mesh)
        beta_15_g_s = np.abs((U[1][:] + V[1][:])/mesh)
        beta_20_g_s = np.abs((U[2][:] + V[2][:])/mesh)

        # ========================================================================

        # Plotting results
        # ----------------

        # Path for saving the figure
        plot_path = os.path.join(glob_path, "Results", f"Beta_Pc={press}_Pw={Pow}.pdf")

        # Mesh normalization
        normalized_mesh = mesh/R

        # Plotting beta
        plt.figure()
        plt.plot(normalized_mesh,beta_10_g_s)
        plt.plot(normalized_mesh,beta_15_g_s)
        plt.plot(normalized_mesh,beta_20_g_s)
        plt.legend(["mdot= 10 g/s", "mdot= 15 g/s", "mdot= 20 g/s"])
        plt.xlabel("x/R [m]")
        plt.ylabel("beta [m/s]")
        plt.title(f"Pc = {press} [mbar] | Pw = {Pow} [kW]")
        plt.grid(True)  
        plt.savefig(plot_path, format='pdf', dpi=600, bbox_inches='tight', transparent=True)

        #plt.figure()
        #plt.plot(normalized_mesh,np.abs(U[0][:]))
        #plt.show()

        #plt.figure()
        #plt.plot(normalized_mesh,np.abs(V[0][:]))
        #plt.show()












