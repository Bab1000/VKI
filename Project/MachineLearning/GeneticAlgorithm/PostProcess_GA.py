import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Enables 3D plotting
import os 

# ===========================================================================================

# -------------------
# | READING RESULTS |
# -------------------

# Path to CSV file
Results_path = "/home/jpe/VKI/Project/MachineLearning/GeneticAlgorithm/GA_results_exp_campaign.csv"

# Reading data
df = pd.read_csv(Results_path, sep=" ")

#print("mdot uniques :", df["mdot"].unique())
#print("Pc uniques :", df["Pc"].unique())

# ========================================================================
# ---------------------
# | TARGETS FOR PLOTS | 
# ---------------------

target_Pc = [1500, 5000, 10000]
target_mdot = [10, 16, 20]

# ========================================================================
# ----------------------- 
# | PLOTTING q_tot vs T | 
# -----------------------

for stat_pressure in target_Pc:

    df_sorted_pressure = df[df["Pc"] == stat_pressure]

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    mdot_colors = {10: 'tab:blue', 16: 'tab:orange', 20: 'tab:green'}  
    line_styles = "--"  
    line_width = 2  

    for i, mdot_i in enumerate(target_mdot):

        df_sorted_mdot = df_sorted_pressure[df_sorted_pressure["mdot"] == mdot_i]
        ax = axes[i]

        ax.scatter(df_sorted_mdot["Tin"], df_sorted_mdot["Q_exp"], 
                   color=mdot_colors[mdot_i], marker='D', s=80,
                   edgecolor='black', linewidths=2, alpha=0.9, zorder=3,  
                   label="Experimental Data")

        ax.plot(df_sorted_mdot["Tin"], df_sorted_mdot["Q_pred"], 
                marker = 'o',
                linestyle=line_styles, color=mdot_colors[mdot_i],
                linewidth=line_width, label="Predicted Data")

        ax.set_ylabel("q_tot (kW/m²)", fontsize=12)
        ax.set_title(f"mdot = {mdot_i} g/s", fontsize=12)
        ax.grid(True, linestyle='--', linewidth=0.7)
        ax.legend(fontsize=10, loc="upper left", frameon=True)

    axes[-1].set_xlabel("Temperature T (K)", fontsize=12)

    plot_path = os.path.join("/home/jpe/VKI/Project/MachineLearning/GeneticAlgorithm/Results", f"HF_Pc={stat_pressure}Pa.jpeg")
    plt.savefig(plot_path, format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
    plt.close()

# ========================================================================
# -----------------------------
# | PLOTTING GammaN VS GammaO | 
# -----------------------------

for stat_pressure in target_Pc:

    df_sorted_pressure = df[df["Pc"] == stat_pressure]

    # Create 3 subplots in 1 row, each 3D
    fig = plt.figure(figsize=(18, 6))
    mdot_colors = {10: 'tab:blue', 16: 'tab:orange', 20: 'tab:green'}

    for i, mdot_i in enumerate(target_mdot):
        df_sorted_mdot = df_sorted_pressure[df_sorted_pressure["mdot"] == mdot_i]

        # Add a 3D subplot
        ax = fig.add_subplot(1, 3, i+1, projection='3d')

        # Scatter plot
        ax.scatter(
            df_sorted_mdot["gammaN"],
            df_sorted_mdot["gammaO"],
            df_sorted_mdot["Tin"],
            color=mdot_colors[mdot_i],
            marker='o',
            s=60,
            edgecolor='black',
            linewidths=1,
            alpha=0.8
        )

        # Axis labels and title
        ax.set_xlabel("GammaN", fontsize=10)
        ax.set_ylabel("GammaO", fontsize=10)
        ax.set_zlabel("Tin (K)", fontsize=10)
        ax.set_title(f"mdot = {mdot_i} g/s", fontsize=12)

        # Optional: Adjust the view angle
        ax.view_init(elev=30, azim=135)

    # Main title
    fig.suptitle(f"3D Gamma vs Tin | Pc = {stat_pressure} mbar", fontsize=16)

    # Save the figure
    plot_path = os.path.join(
        "/home/jpe/VKI/Project/MachineLearning/GeneticAlgorithm/Results",
        f"3x3D_Gamma_Pc={stat_pressure}Pa.jpeg"
    )
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)  # leave space for suptitle
    plt.savefig(plot_path, format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
    plt.close()