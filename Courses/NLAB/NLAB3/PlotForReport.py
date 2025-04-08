import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re

# --------------
# | QUESTION 2 |
# --------------

# Gather results for plots in mesh study path and for mesh = 10023

# Read the data
file_path = "/home/jpe/VKI/Courses/NLAB/NLAB3/MeshStudy/10023/history.csv"  # Change this to your file path
df = pd.read_csv(file_path, skipinitialspace=True)
df.columns = df.columns.str.strip()

resP = df["rms[P]"]
resU = df["rms[U]"]
resV = df["rms[V]"]
CFL = df["Avg CFL"]
iter = df["Inner_Iter"]
 
# Plotting results
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(iter, resP, label="Pressure",linewidth=2,)
ax1.plot(iter, resU, label="U",linewidth=2,)
ax1.plot(iter, resV, label="V",linewidth=2,)
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Residuals')
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(iter, CFL, label="CFL", color="red",linewidth=2, linestyle="dashed")
ax2.set_ylabel('CFL')

fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.title('Residuals and CFL vs Iterations')
plt.savefig("/home/jpe/VKI/Courses/NLAB/NLAB3/Figures/Res_CFL_vs_iter.png", format='png', dpi=300, bbox_inches='tight')
 
# =========================================================================================

# --------------
# | QUESTION 5 |
# --------------

# Mesh study

# Numbers of elements in the meshes studied
nbr_elem = [3576, 7062, 10023]

# Creating storing vectors
Cl_val_vec = []
Cd_val_vec = []

# Loop on the number of mesh studied
for mesh_nbr in nbr_elem:

    # Path of the SU2 simulation
    Sim_path = f"/home/jpe/VKI/Courses/NLAB/NLAB3/MeshStudy/{mesh_nbr}/"

    # Forces file
    Force_file_path = Sim_path + "forces_breakdown.dat"

    # Initialisation des valeurs
    total_Cl = None
    total_Cd = None

    with open(Force_file_path, "r", encoding="utf-8") as file:
        for line in file:
            match_Cl = re.search(r"Total CL:\s*([-+]?\d*\.\d+)", line)
            match_Cd = re.search(r"Total CD:\s*([-+]?\d*\.\d+)", line)

            if match_Cl:
                total_Cl = float(match_Cl.group(1))

            if match_Cd:
                total_Cd = float(match_Cd.group(1))

            if total_Cl is not None and total_Cd is not None:
                break

    if total_Cl is not None and total_Cd is not None:
        print(f"For element number = {mesh_nbr}: CL = {total_Cl} | CD = {total_Cd}")
    else:
        print(f"No CL values found in {Force_file_path}")

    # Ajout aux listes
    Cl_val_vec.append(total_Cl)
    Cd_val_vec.append(total_Cd)

# Plotting results
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(nbr_elem, Cl_val_vec,linewidth=2, label="CL values")
ax1.set_xlabel("Number of elements in the mesh")
ax1.set_ylabel("CL Value")
ax1.tick_params(axis='y')
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(nbr_elem, Cd_val_vec, linewidth=2, color="orange", label="CD values")
ax2.set_ylabel("CD Value")
ax2.tick_params(axis='y')

fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.title('Aerodynamic coefficients vs number of elements in the mesh')
plt.xticks(nbr_elem, ["8116", "12178", "15523"])
plt.savefig("/home/jpe/VKI/Courses/NLAB/NLAB3/Figures/Mesh_study.png", format='png', dpi=300, bbox_inches='tight')


# =========================================================================================

# --------------
# | QUESTION 6 |
# --------------

# Path of the SU2 simulation for CFL strategy study
Sim_path = f"/home/jpe/VKI/Courses/NLAB/NLAB3/CFLStrategyStudy/"

history_file_path = Sim_path + "history_adapt.csv"

df = pd.read_csv(history_file_path, skipinitialspace=True)
df.columns = df.columns.str.strip()

resP = df["rms[P]"]
resU = df["rms[U]"]
resV = df["rms[V]"]
CFL = df["Avg CFL"]
iter = df["Inner_Iter"]
 
# Plotting results
# Create figure and axes
fig, ax1 = plt.subplots(figsize=(10, 6))

ax2 = ax1.twinx()
ax2.plot(iter, CFL, label="CFL", color="red", linestyle="dashed", linewidth=2, alpha=0.3, zorder=1)
ax2.set_ylabel('CFL')

# Plot residuals (foreground)
ax1.plot(iter, resP, label="Pressure", linewidth=2, zorder=2)
ax1.plot(iter, resU, label="U", linewidth=2, zorder=2)
ax1.plot(iter, resV, label="V", linewidth=2, zorder=2)
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Residuals')
ax1.grid(True)

# Add legend
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

# Title and save the plot
plt.title('Residuals and CFL vs Iterations')
plt.savefig("/home/jpe/VKI/Courses/NLAB/NLAB3/Figures/Adaptive_CFL.png", format='png', dpi=300, bbox_inches='tight')


# =========================================================================================

# --------------
# | QUESTION 7 |
# --------------


# Path of the SU2 simulation for CFL strategy study
Sim_path = f"/home/jpe/VKI/Courses/NLAB/NLAB3/UnsteadyCorentin/"

history_file_path = Sim_path + "history.csv"

df = pd.read_csv(history_file_path, skipinitialspace=True)
df.columns = df.columns.str.strip()

CD = df["CD"]
CL = df["CL"]
iter = df["Time_Iter"]

plt.figure(figsize=(10,6))
plt.plot(iter[30:]/20,CD[30:])
plt.xlabel("Time [s]")
plt.ylabel("CD [-]")
plt.grid(True)
plt.savefig("/home/jpe/VKI/Courses/NLAB/NLAB3/Figures/CD_unsteady.png", format='png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10,6))
plt.plot(iter[30:]/20,CL[30:])
plt.xlabel("Time [s]")
plt.ylabel("CL [-]")
plt.grid(True)
plt.savefig("/home/jpe/VKI/Courses/NLAB/NLAB3/Figures/CL_unsteady.png", format='png', dpi=300, bbox_inches='tight')


