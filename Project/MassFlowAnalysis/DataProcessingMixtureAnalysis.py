import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Utils_DataProcessingMixtureAnalysis import *

# Target pressure and power for mass flow analysis
targ_p = [15,50,100]
targ_P = [150, 200, 250, 300, 350]

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
            pres_pow_pattern_name = f"Pc={press}_Pw={Pow}"
            massflow_pattern_name = "mdot=*"
            surface_pattern_name = "*surface.dat"

            pres_pow_folder_path = os.path.join(g_path, pres_pow_pattern_name)
            surface_general_path = os.path.join(pres_pow_folder_path, massflow_pattern_name, "**", surface_pattern_name)

            # Locate the flowfield data files
            list_path = glob.glob(surface_general_path, recursive=True)

            for file in list_path:
                mdot = extract_mdot(file)
                temp = extract_temperature(file)
                
                if mdot is not None and temp is not None:
                    q_tot = extract_q_tot(file)  # Extract q_tot
                    all_results.append((mdot, temp, q_tot, label))

    # Convert results into a DataFrame for processing
    df = pd.DataFrame(all_results, columns=["mdot", "Temperature", "q_tot", "Label"])

    # Apply rounding to mdot values
    df["mdot"] = df["mdot"].apply(round_mdot)

    # Sort values by Label, mdot, and Temperature
    df = df.sort_values(by=["Label", "mdot", "Temperature"])

    # ========================================================================
    # ------------------------------
    # | GENERIC DATA LOADING LOGIC |
    # ------------------------------

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
    # ------------ 
    # | PLOTTING | 
    # ------------ 

    # Initialize plot with a larger figure size
    plt.figure(figsize=(12, 7))

    # Define line styles for each air type
    line_styles = {"air_5": "-", "air_7": "--"}  # Solid for air_5, dashed for air_7

    # Define colors for each mdot value (consistent across both datasets)
    mdot_colors = {10: 'tab:blue', 16: 'tab:orange', 20: 'tab:green'}  

    # Increase line width for simulation data
    line_width = 2

    # First, plot simulation data (so that experimental data will be on top)
    for label in df["Label"].unique():
        df_label = df[df["Label"] == label]

        for mdot in [10, 16, 20]:  # Loop over mdot values
            df_mdot = df_label[df_label["mdot"] == mdot]

            if not df_mdot.empty:
                # Plot simulation data first
                plt.plot(df_mdot["Temperature"], np.abs(df_mdot["q_tot"]) / 1000, 
                        linestyle=line_styles[label], 
                        marker='o', markersize=7,  
                        color=mdot_colors[mdot], 
                        linewidth=line_width,  
                        label=f"{label} - mdot {mdot} g/s")

    # Now, plot experimental data **last** so it's on top
    for mdot in [10, 16, 20]:
        if len(T_mdot_data[mdot]) > 0 and len(exp_HF_data[mdot]) > 0:
            plt.scatter(T_mdot_data[mdot], np.abs(exp_HF_data[mdot]), 
                        color=mdot_colors[mdot], marker='D', 
                        s=80,  # Increased marker size for better visibility
                        edgecolor='black', linewidths=2,  # Stronger border for contrast
                        alpha=0.9,  # Full opacity for clear distinction
                        zorder=3,  # Ensures experimental points are drawn **above** all others
                        label=f"Exp. val. - mdot {mdot} g/s")

    # Configure labels, title, and grid
    plt.xlabel("Temperature T (K)", fontsize=14)
    plt.ylabel("Total heat flux q_tot (kW/m²)", fontsize=14)
    plt.title(f"Heat flux vs Temperature Pc={press}mbar (Stagline vs Experiment)", fontsize=16)
    plt.legend(fontsize=12, loc="upper left", frameon=True)  
    plt.grid(True, linestyle='--', linewidth=0.7)
    comparison_plot_path = os.path.join("/home/jpe/VKI/Project/MassFlowAnalysis/Results_HFComp", f"HFComparison_Pc={press}.pdf")
    plt.savefig(comparison_plot_path, format='pdf', dpi=600, bbox_inches='tight', transparent=True)

    all_results.clear()

