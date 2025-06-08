import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from Utils_PostProcessingBeta import extract_BL_edge,FlowfieldReader,extract_pressure,extract_h_wall,extract_h_edge,extract_R
from scipy.interpolate import interp1d

case = "CASE_9"
M1 = 6

flowfield_file_path = f"/home/jpe/VKI/Courses/CSA/Stagline_simulations_Vitto/{case}/Mach/M1_{M1}/{case}_Air5_SMB_MPP_2R/Air5_SMB_MPP_2R_flowfield.dat"
surface_file_path = f"/home/jpe/VKI/Courses/CSA/Stagline_simulations_Vitto/{case}/Mach/M1_{M1}/{case}_Air5_SMB_MPP_2R/Air5_SMB_MPP_2R_surface.dat"
mesh_file_path = f"/home/jpe/VKI/Courses/CSA/Stagline_simulations_Vitto/{case}/Mach/M1_{M1}/{case}_Air5_SMB_MPP_2R/mesh.dat"

# Extracting BL_thickness from the surface file
BL_edge = extract_BL_edge(surface_file_path)

# Extracting flowfield data
mesh,U,V,T,pc = FlowfieldReader(flowfield_file_path)

Pstag = extract_pressure(surface_file_path)

h_wall = extract_h_wall(surface_file_path)

h_edge = extract_h_edge(surface_file_path)

R = extract_R(mesh_file_path)

# ========================================================================
 
# ----------------------------------------------
# | BETA PROFILE COMPUTATION + BETA AT BL EDGE |
# ----------------------------------------------


plot_path_beta = f"/home/jpe/VKI/Courses/CSA/Plots_Vitto/{case}/Mach/M1_{M1}/{case}_Air5_SMB_MPP_2R"
os.makedirs(plot_path_beta, exist_ok=True)


# Compute velocity gradient profile
beta_profile = np.abs((U + V) / mesh)

# BL edge position
BL_edge_pos = R + BL_edge

# Compute velocity gradient at BL edge
beta_interp = interp1d(mesh, beta_profile)
beta_edge = beta_interp(BL_edge_pos)

# Normilized mesh
normalized_mesh = mesh/R

plt.figure(figsize=(8, 6))
plt.plot(normalized_mesh, beta_profile, label=f"Velocity gradient profile")
plt.scatter(BL_edge_pos/R, beta_edge, s=100,color="red", label=f"Beta BL edge : {beta_edge:.3f} [1/s]")
plt.xlabel("x/R [-]", fontsize=18)
plt.ylabel("beta [1/s]", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True)
plt.savefig(plot_path_beta + "/beta.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
plt.close()

print(f"Stagnation pressure : {Pstag:.3f} [Pa]")
print(f"Wall enthalpy : {h_wall:.3f} [Pa]")
print(f"BL edge enthalpy : {h_edge:.3f} [Pa]")
print(f"BL edge position : {BL_edge_pos:.5f} [m]")
print(f"Beta at BL edge : {beta_edge:.2f} [1/s]")

# ========================================================================
 
# -----------------------------------
# | TEMPERATURE PROFILE COMPUTATION |
# -----------------------------------

# Normilized mesh
normalized_mesh = mesh/R

# find max temperature
max_T = np.max(T)

# find index on the mesh
idx = np.argmax(T)

# Compute shoch width from the probe surface 
pos_shock = mesh[idx]
width_shock = pos_shock - R

print(f"Shock width : {width_shock:.5f} [m]")


plt.figure(figsize=(8, 6))
plt.plot(normalized_mesh, T, label="Temperature profile")
plt.xlabel("x/R [-]", fontsize=18)  # Taille augmentée ici
plt.ylabel("T [K]", fontsize=18)    # Taille augmentée ici
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True)
plt.savefig(plot_path_beta + "/Temperature_profile.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
plt.close()
