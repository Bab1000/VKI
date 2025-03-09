import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Utils_DataProcessingRebuildingAnalysis import *

# Target pressure and power for mass flow analysis
targ_p = [15,50,100]
targ_P = [150, 200, 250, 300, 350]

# Probe radius
R = 0.025

# Paths for data processing
glob_path = {
    "air_5": "/home/jpe/VKI/Project/MassFlowAnalysis/air_5_sim",
    "air_7": "/home/jpe/VKI/Project/MassFlowAnalysis/air_7_sim"
}

# Store results
all_results = []

# Loop over pressure and power values
for press in targ_p:
    for Pow in targ_P:
        for label, g_path in glob_path.items():

    # ========================================================================
    # -------------------------------------------
    # | GATHERING HF/H/BL EDGE FROM SYRFACE FILE|
    # -------------------------------------------

            # Building pattern for folders and files
            press_pow_pattern_name = f"Pc={press}_Pw={Pow}"
            massflow_pattern_name = "mdot=*"
            surface_pattern_name = "*surface.dat"
            flowfield_pattern_name = "*flowfield.dat"
            
            # Building the path for data folder
            press_pow_folder_path = os.path.join(g_path, press_pow_pattern_name)
            # Building the path for surface.dat file
            surface_general_path = os.path.join(press_pow_folder_path, massflow_pattern_name, "**", surface_pattern_name)
            # Building the path for surface.dat file
            flowfield_general_path = os.path.join(press_pow_folder_path, massflow_pattern_name, "**", flowfield_pattern_name)

            # Locate the surface data files
            list_path_surface = glob.glob(surface_general_path, recursive=True)

            # Locate the flowfield data files
            list_path_flowfield = glob.glob(flowfield_general_path, recursive=True)

            counter = 0
            for file in list_path_surface:

                # Extracting simulation massflow/temperature from path name
                mdot = extract_mdot(file)
                temp = extract_temperature(file)

                # Extracting q_tot, h_edge and BL_thickness from the surface file
                q_tot = extract_q_tot(file)  
                h_edge = extract_h_edge(file)
                BL_edge = extract_BL_edge(file)
                u_edge = extract_u_edge(file)
                density_edge = extract_density_edge(file)
                heat_transfer_coeff = extract_heat_transfer_coeff(file)
                
                # building flowflied file path
                flowfield_file_path = list_path_flowfield[counter]

                # Extracting flowfield data
                mesh,U,V,T,pc = FlowfieldReader(flowfield_file_path)

                # Position of the BL edge
                BL_edge_pos = R + BL_edge

                # Find the indices of the closest points surrounding the BL thickness
                idx_below = np.where(mesh < BL_edge_pos)[0][-1]  # Last point below BL
                idx_above = np.where(mesh > BL_edge_pos)[0][0]    # First point above BL
                
                # Average velocities
                avg_U = np.abs(U[idx_above] + U[idx_below])/2 
                avg_V = np.abs(V[idx_above] + V[idx_below])/2

                # Average static pressure
                avg_pc = np.abs(pc[idx_above] + pc[idx_below])/2

                # Delta T at BL edge
                #delta_T = np.abs(T[idx_above] - T[idx_below])

                # BL edge heat flux 
                #q_tot = heat_transfer_coeff * delta_T

                # BL edge velocity gradient computation
                beta = (avg_U + avg_V) / BL_edge_pos

                # BL edge total pressure computation
                p0 = avg_pc + 0.5 * density_edge * (avg_U**2 + avg_V**2)

                if mdot is not None and temp is not None and q_tot is not None and h_edge is not None and beta is not None:
                    all_results.append((mdot, temp, q_tot, h_edge, beta, p0, label))

                counter = counter + 1


    # ========================================================================
    # ---------------------- 
    # | BUILDING DATAFRAME | 
    # ----------------------

    # Convert results into a DataFrame for processing
    df = pd.DataFrame(all_results, columns=["mdot", "Temperature", "q_tot", "h_edge", "beta", "p0", "Label"])

    # Apply rounding to mdot values
    df["mdot"] = df["mdot"].apply(round_mdot)

    # Sort values by Label, mdot, and Temperature
    df = df.sort_values(by=["Label", "mdot", "Temperature"])

    # ========================================================================
    # -----------------------------------
    # | EXPERIMENTAL DATA LOADING LOGIC |
    # -----------------------------------

    # Path to the CSV
    CSV_path = "/home/jpe/VKI/Project/MassFlowAnalysis/tests_Justin.xlsx"

    # Define the columns to check for NaN values
    columns_to_check = ["Pressure[mbar]", "massflow [g/s]", "Power[kW]", "HeatFlux(HS50mm)[kW/m2]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)"]

    # Load data from CSV
    pressure, massflow, power, heat_flux, pitot, temperature = CSVReader(CSV_path, columns_to_check)

    # Define target mdot values
    target_mdot_values = [10, 16, 20]

    # Define tolerance range for filtering
    range_tolerance = 3

    # Dictionary to store experimental heat flux and temperature for each mdot
    exp_HF_data = {mdot: [] for mdot in target_mdot_values}
    T_mdot_data = {mdot: [] for mdot in target_mdot_values}

    # Filter experimental data based on target pressure and mdot values
    for i, exp_press in enumerate(pressure):
        if press - range_tolerance <= exp_press <= press + range_tolerance:
            for mdot in target_mdot_values:
                if mdot - range_tolerance <= massflow[i] <= mdot + range_tolerance:
                    exp_HF_data[mdot].append(heat_flux[i])
                    T_mdot_data[mdot].append(temperature[i])

    # Convert lists to NumPy arrays for better handling
    for mdot in target_mdot_values:
        exp_HF_data[mdot] = np.array(exp_HF_data[mdot])
        T_mdot_data[mdot] = np.array(T_mdot_data[mdot])



    # ========================================================================
    # ------------------ 
    # | PLOTTING q_tot | 
    # ------------------

    # Initialize subplots (3 rows, 1 column, sharing X-axis)
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)

    # Define line styles for each air type
    line_styles = {"air_5": "-", "air_7": "--"}  # Solid for air_5, dashed for air_7

    # Define colors for each mdot value (consistent across both datasets)
    mdot_colors = {10: 'tab:blue', 16: 'tab:orange', 20: 'tab:green'}  

    line_width = 2  # Increase line width for simulation data

    # Loop through each mdot value and plot in its corresponding subplot
    for i, mdot in enumerate([10, 16, 20]):
        ax = axes[i]  # Select the corresponding subplot
        
        # Plot simulation data first (so experimental data appears on top)
        for label in df["Label"].unique():
            df_label = df[df["Label"] == label]
            df_mdot = df_label[df_label["mdot"] == mdot]
            
            if not df_mdot.empty:
                ax.plot(df_mdot["Temperature"], np.abs(df_mdot["q_tot"]) / 1000, 
                        linestyle=line_styles[label], 
                        marker='o', markersize=7,  
                        color=mdot_colors[mdot], 
                        linewidth=line_width,  
                        label=f"{label} - Simulation")
        
        # Plot experimental data (scatter plot, on top of simulation curves)
        if len(T_mdot_data[mdot]) > 0 and len(exp_HF_data[mdot]) > 0:
            ax.scatter(T_mdot_data[mdot], np.abs(exp_HF_data[mdot]), 
                    color=mdot_colors[mdot], marker='D', 
                    s=80, edgecolor='black', linewidths=2,  
                    alpha=0.9, zorder=3,  
                    label="Experimental Data")
        
        # Configure subplot
        ax.set_ylabel("q_tot (kW/m²)", fontsize=12)
        ax.set_title(f"mdot = {mdot} g/s", fontsize=12)
        ax.grid(True, linestyle='--', linewidth=0.7)
        ax.legend(fontsize=10, loc="upper left", frameon=True)

    # Global X-label
    axes[-1].set_xlabel("Temperature T (K)", fontsize=12)

    # Save the figure
    comparison_plot_path = os.path.join("/home/jpe/VKI/Project/MassFlowAnalysis/Results_Rebuilding", f"HFComparison_Pc={press}.jpeg")
    plt.savefig(comparison_plot_path, format='jpeg', dpi=600, bbox_inches='tight', transparent=True)



    # ========================================================================
    # -------------------
    # | PLOTTING h_edge | 
    # -------------------

    fig_hedge, axes_hedge = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    
    for i, mdot in enumerate([10, 16, 20]):
        ax = axes_hedge[i]
        
        for label in df["Label"].unique():
            df_label = df[df["Label"] == label]
            df_mdot = df_label[df_label["mdot"] == mdot]
            
            if not df_mdot.empty:
                ax.plot(df_mdot["Temperature"], np.abs(df_mdot["h_edge"]) / 1e6, linestyle='-', marker='s', markersize=7, label=f"{label} - Simulation")
                
        ax.set_ylabel("h_edge (MJ/kg)")
        ax.set_title(f"mdot = {mdot} g/s")
        ax.grid(True, linestyle='--', linewidth=0.7)
        ax.legend()
    
    axes_hedge[-1].set_xlabel("Temperature T (K)")
    plt.savefig(f"/home/jpe/VKI/Project/MassFlowAnalysis/Results_Rebuilding/h_edgeComparison_Pc={press}.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)

    # ========================================================================
    # -------------------
    # | PLOTTING beta | 
    # -------------------

    fig_beta, axes_beta = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    
    for i, mdot in enumerate([10, 16, 20]):
        ax = axes_beta[i]
        
        for label in df["Label"].unique():
            df_label = df[df["Label"] == label]
            df_mdot = df_label[df_label["mdot"] == mdot]
            
            if not df_mdot.empty:
                ax.plot(df_mdot["Temperature"], np.abs(df_mdot["beta"]), linestyle='-', marker='s', markersize=7, label=f"{label} - Simulation")
                
        ax.set_ylabel("beta (1/s)")
        ax.set_title(f"mdot = {mdot} g/s")
        ax.grid(True, linestyle='--', linewidth=0.7)
        ax.legend()
    
    axes_beta[-1].set_xlabel("Temperature T (K)")
    plt.savefig(f"/home/jpe/VKI/Project/MassFlowAnalysis/Results_Rebuilding/betaComparison_Pc={press}.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)

    # ========================================================================
    # -------------------
    # | PLOTTING p0 | 
    # -------------------

    fig_p0, axes_p0 = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    
    for i, mdot in enumerate([10, 16, 20]):
        ax = axes_p0[i]
        
        for label in df["Label"].unique():
            df_label = df[df["Label"] == label]
            df_mdot = df_label[df_label["mdot"] == mdot]
            
            if not df_mdot.empty:
                ax.plot(df_mdot["Temperature"], np.abs(df_mdot["p0"]), linestyle='-', marker='s', markersize=7, label=f"{label} - Simulation")
                
        ax.set_ylabel("p0 (Pa)")
        ax.set_title(f"mdot = {mdot} g/s")
        ax.grid(True, linestyle='--', linewidth=0.7)
        ax.legend()
    
    axes_p0[-1].set_xlabel("Temperature T (K)")
    plt.savefig(f"/home/jpe/VKI/Project/MassFlowAnalysis/Results_Rebuilding/p0Comparison_Pc={press}.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
    
    all_results.clear()
