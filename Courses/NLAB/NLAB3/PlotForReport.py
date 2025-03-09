import numpy as np
import matplotlib.pyplot as plt
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')  
plt.rc('axes', unicode_minus=False)

#%% Inlet Profile (Reference Mesh)
inlet = np.loadtxt('inlet.dat', skiprows=5)

plt.plot(inlet[:,1], inlet[:,3], linewidth=1.5)
plt.xlabel('y [m]', fontsize=18)
plt.ylabel('U(0,y) [m/s]', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True)
plt.show()

#%% Convegerce Reference Mesh (hi/h1=2): Residuals vs Iterations vs CFL
history = np.loadtxt('history.txt', skiprows=2)

fig, ax1 = plt.subplots()
# Plot Residuals
l1 = ax1.plot(history[:,0], history[:,1], linewidth=1.5, label="Pressure")
l2 = ax1.plot(history[:,0], history[:,2], linewidth=1.5, label="U-component")
l3 = ax1.plot(history[:,0], history[:,3], linewidth=1.5, label="V-component")
# Labels and styling for Residuals
ax1.set_xlabel('Iterations', fontsize=18)
ax1.set_ylabel('Residuals', fontsize=18)
ax1.tick_params(axis='both', labelsize=18)
ax1.grid(True)
# Plot CFL
ax2 = ax1.twinx()
l4 = ax2.axhline(y=100, color='r', linestyle='-', linewidth=1.5, label='CFL')
# Labels for CFL
ax2.set_ylabel('CFL', fontsize=18)
ax2.tick_params(axis='y', labelsize=18)
# Combine legends from both axes
lines = l1 + l2 + l3 + [l4] 
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, fontsize=18)
plt.show()

#%% Grid Convergence Study
CL = np.array([0.01364198587, 0.01007849909, 0.01043424399]) # Splitting refinement
CD = np.array([5.657852877, 5.617917356, 5.585568726])       # Splitting refinement
# CL = np.array([0.01364198587, 0.01007849909, 0.01017546143])  # Doubling Ncells
# CD = np.array([5.657852877, 5.617917356, 5.621656999])
x_values = [4, 2, 1]

# Plot CL
fig, ax = plt.subplots()
ax.semilogy(x_values, CL, 'o-', linewidth=1.5)
ax.set_xlabel('$h_{i}/h_{1}$', fontsize=18)
ax.set_ylabel('$C_L$', fontsize=18)
ax.set_xticks(x_values)
ax.set_xticklabels(['4', '2', '1'], fontsize=16)
ax.set_yticks(CL)
ax.set_yticklabels([f'{v:.4f}' for v in CL], fontsize=18)
ax.yaxis.set_minor_locator(plt.NullLocator())  
ax.grid(True)
plt.show()

# Plot CD
fig, ax = plt.subplots()
ax.semilogy(x_values, CD, 'o-', linewidth=1.5)
ax.set_xlabel('$h_{i}/h_{1}$', fontsize=18)
ax.set_ylabel('$C_D$', fontsize=18)
ax.set_xticks(x_values)
ax.set_xticklabels(['4', '2', '1'], fontsize=18)
ax.set_yticks(CD)
ax.set_yticklabels([f'{v:.4f}' for v in CD], fontsize=18)
ax.yaxis.set_minor_locator(plt.NullLocator())  
ax.grid(True)
plt.show()

#%% Adaptive CFL - h2 (Reference Mesh)
path = r"C:\Users\User\Desktop\VKI - RM\Lectures\NLAB1\HM_3\CFL_Adaptive\h2"
os.chdir(path)

history_CFL_adpt = np.loadtxt('history.txt', skiprows=2)

fig, ax1 = plt.subplots()
# Plot Residuals
l1 = ax1.plot(history_CFL_adpt[:,0], history_CFL_adpt[:,1], linewidth=1.5, label="Pressure")
l2 = ax1.plot(history_CFL_adpt[:,0], history_CFL_adpt[:,2], linewidth=1.5, label="U-component")
l3 = ax1.plot(history_CFL_adpt[:,0], history_CFL_adpt[:,3], linewidth=1.5, label="V-component")
# Labels and styling for Residuals
ax1.set_xlabel('Iterations', fontsize=18)
ax1.set_ylabel('Residuals', fontsize=18)
ax1.tick_params(axis='both', labelsize=18)
ax1.grid(True)
# Plot CFL
ax2 = ax1.twinx()
l4 = ax2.plot(history_CFL_adpt[:,0], history_CFL_adpt[:,4], linewidth=1.5, alpha=0.3, label="CFL")
# Labels for CFL
ax2.set_ylabel('CFL', fontsize=18)
ax2.tick_params(axis='y', labelsize=18)
# Combine legends from both axes
lines = l1 + l2 + l3 + l4 
labels = [line.get_label() for line in lines]
legend = ax1.legend(lines, labels, fontsize=18, loc='best', frameon=True)
legend.set_zorder(100)  
plt.show()

#%% Unsteady 


