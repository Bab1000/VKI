import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Target pressure and power for the mass flow analysis
targ_p = [50, 100]
targ_P = [150, 200, 250, 300, 350]

# Paths for data processing
glob_path = {
    "air_5": "/home/jpe/VKI/Project/MassFlowAnalysis/air_5_sim",
    "air_7": "/home/jpe/VKI/Project/MassFlowAnalysis/air_7_sim"
}

# Probe radius [m]
R = 0.025

# Couleurs uniformes pour la comparaison
colors = ["tab:blue", "tab:orange", "tab:green"]  # Bleu, orange, vert

# Loop over pressure and power
for press in targ_p:
    for Pow in targ_P:
        results = {}

        for label, g_path in glob_path.items():
            pres_pow_pattern_name = f"Pc={press}_Pw={Pow}"
            massflow_pattern_name = "mdot=*"
            flowfield_pattern_name = "*flowfield.dat"

            pres_pow_folder_path = os.path.join(g_path, pres_pow_pattern_name)
            flowfiled_general_path = os.path.join(pres_pow_folder_path, massflow_pattern_name, "**", flowfield_pattern_name)

            # Locate the flowfield data files
            list_path = glob.glob(flowfiled_general_path, recursive=True)

            # Sorting the files by massflow values
            def extract_mdot_value(filepath):
                match = re.search(r"mdot=([\d\.]+)", filepath)
                return float(match.group(1)) if match else float('inf')

            list_path = sorted(list_path, key=extract_mdot_value)

            # Storing matrices
            U, V, mesh = [], [], None

            for flowfield_path in list_path:
                df = pd.read_csv(flowfield_path, delim_whitespace=True, header=None, skiprows=4)
                u = df.iloc[:, 1].tolist()
                v = df.iloc[:, 2].tolist()
                mesh = df.iloc[:, 0].tolist()

                U.append(u)
                V.append(v)

            # Transform to arrays
            U, V, mesh = np.array(U), np.array(V), np.array(mesh)

            # Compute velocity gradient
            beta_10_g_s = np.abs((U[0][:] + V[0][:]) / mesh)
            beta_15_g_s = np.abs((U[1][:] + V[1][:]) / mesh)
            beta_20_g_s = np.abs((U[2][:] + V[2][:]) / mesh)

            # Store results
            results[label] = {
                "mesh": mesh / R,  # Normalized mesh
                "beta": [beta_10_g_s, beta_15_g_s, beta_20_g_s]  # List to keep all beta values
            }

            # ========================================================================
            # Plot individual results (like before)
            plot_path = os.path.join(g_path, "Results", f"Beta_Pc={press}_Pw={Pow}.pdf")

            plt.figure(figsize=(8, 6))
            for beta, color, mdot in zip(results[label]["beta"], colors, ["10 g/s", "15 g/s", "20 g/s"]):
                plt.plot(results[label]["mesh"], beta, color=color, label=f"mdot={mdot}")

            plt.xlabel("x/R [-]")
            plt.ylabel("beta [1/s]")
            plt.title(f"{label} | Pc = {press} [mbar] | Pw = {Pow} [kW]")
            plt.legend()
            plt.grid(True)
            plt.savefig(plot_path, format='pdf', dpi=600, bbox_inches='tight', transparent=True)
            plt.close()

        # ========================================================================
        # Plot comparative results (air_5 vs air_7)
        plt.figure(figsize=(8, 6))

        for i, (mdot, color) in enumerate(zip(["10 g/s", "16 g/s", "20 g/s"], colors)):
            plt.plot(results["air_5"]["mesh"], results["air_5"]["beta"][i], linestyle="--", color=color, label=f"air_5 mdot={mdot}")
            plt.plot(results["air_7"]["mesh"], results["air_7"]["beta"][i], linestyle="-", color=color, label=f"air_7 mdot={mdot}")

        plt.xlabel("x/R [-]")
        plt.ylabel("beta [1/s]")
        plt.title(f"Beta comparison | Pc = {press} [mbar] | Pw = {Pow} [kW]")
        plt.legend()
        plt.grid(True)

        # Save the comparative plot
        comparison_plot_path = os.path.join("/home/jpe/VKI/Project/MassFlowAnalysis/Results_MixtComp", f"BetaComparison_Pc={press}_Pw={Pow}.pdf")
        plt.savefig(comparison_plot_path, format='pdf', dpi=600, bbox_inches='tight', transparent=True)
        plt.close()
