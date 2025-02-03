# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 11:52:32 2022

@author: mendez, ratz
"""

# Basics needed to create and plot the data
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration for plots  
fontsize = 12
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['legend.fontsize'] = fontsize-3
plt.rcParams['font.size'] = fontsize

# Define the output folder
Fol_Out = '01_Radial_Basis_Functions'
if not os.path.exists(Fol_Out):
    os.makedirs(Fol_Out)

#%% Define all of the bases and the functions to assemble the basis matrices

# The Gaussian
def Gauss_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r
    phi_r = np.exp(-c_r**2*d**2)
    return phi_r

# The Matern C4 kernel
def Matern_C4_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = np.abs(x-x_r)
    phi_r = np.exp(-c_r*d) * (c_r**2*d**2 + 3*c_r*d + 3)
    return phi_r

# The (custom) C4
def Compact_C4_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r 
    phi_r = (1+d/c_r)**5 * (1-d/c_r)**5
    phi_r[np.abs(d)>c_r] = 0
    return phi_r

#%% Plot the curves

# Generate the sample points and basis locations
x = np.linspace(-1, 1, 1001)
x_b = np.linspace(-1, 1, 5)

# Plot the inverse multiquadratic
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
for i in range(5):
    ax.plot(x, Matern_C4_RBF(x, x_b[i], c_r=10))
ax.set_title('Matern C4 RBF')
ax.set_xlabel('x')
ax.set_ylabel('y')
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'RBF_MaternC4.png', bbox_inches='tight', pad_inches=0)
    

# Plot the Gaussian
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
for i in range(5):
    ax.plot(x, Gauss_RBF(x, x_b[i], c_r=5))
ax.set_title('Gaussian RBF')
ax.set_xlabel('x')
ax.set_ylabel('y')
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'RBF_gaussian.png', bbox_inches='tight', pad_inches=0)


# Plot the C4
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
for i in range(5):
    ax.plot(x, Compact_C4_RBF(x, x_b[i], c_r=0.5))
ax.set_title('C4 RBF')
ax.set_xlabel('x')
ax.set_ylabel('y')
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'RBF_C4.png', bbox_inches='tight', pad_inches=0)

plt.show()
