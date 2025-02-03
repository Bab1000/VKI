# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:15:50 2025

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
Fol_Out = '03_Meshless_PDE'
if not os.path.exists(Fol_Out):
    os.makedirs(Fol_Out)

#%% We seek to solve a boundary value problem in a meshless fashion
# Generate the 'domain'
n_x = 500
n_b = 200
x = np.linspace(0, 1, n_x)
x_b = np.linspace(0, 1, n_b)


#%% Step 1: Create the matrix L, Phi_Gamma, c
# We need Phi, Phi_x and Phi_xx
def C4_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r
    phi_r = (1+d/c_r)**5 * (1-d/c_r)**5
    phi_r[np.abs(d)>c_r] = 0
    return phi_r

def C4_RBF_x(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r
    phi_r=(-10*d*(d-c_r)**4 * (d+c_r)**4) / (c_r**10)
    phi_r[np.abs(d)>c_r] = 0
    return phi_r

def C4_RBF_xx(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r
    phi_r = (-10*(d-c_r)**3 * (d+c_r)**3) * (3*d-c_r) * (3*d+c_r) / (c_r**10)
    phi_r[np.abs(d)>c_r] = 0
    return phi_r


phi = C4_RBF(x, 0.5, c_r=0.05)
phi_x_a = C4_RBF_x(x,0.5,c_r=0.05)
phi_xx_a = C4_RBF_xx(x, 0.5, c_r=0.05)

# # If you want to check the derivatives numerically:
# phi_x_n=np.gradient(phi,x)
# phi_xx_n=np.gradient(phi_x_a,x)

# plt.figure()
# plt.plot(phi_x_a)
# plt.plot(phi_x_n)

# plt.plot(phi_xx_a)
# plt.plot(phi_xx_n)

# plt.show()


# Function that prepares the three derivatives
def PHI_C4_X(x_in, x_b, c_r=0.1):
    n_x=np.size(x_in)
    n_b=len(x_b)
    # Initialize the three bases matrices
    Phi = np.zeros((n_x, n_b))
    Phi_x = np.zeros((n_x, n_b))
    Phi_xx = np.zeros((n_x, n_b))
    
    # Fill the matrices with the C4 and its derivatives
    for j in range(0, n_b):
        Phi[:, j] = C4_RBF(x_in, x_r=x_b[j], c_r=c_r)
        Phi_x[:, j] = C4_RBF_x(x_in, x_r=x_b[j], c_r=c_r)
        Phi_xx[:, j] = C4_RBF_xx(x_in, x_r=x_b[j], c_r=c_r)
 
    return Phi, Phi_x, Phi_xx


Phi, Phi_x, Phi_xx = PHI_C4_X(x, x_b, c_r=0.1)


# Plot the first and second derivative
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
ax.plot(x, Phi_x[:, 150],'b', label='dydx')
ax.plot(x, Phi_xx[:, 100],'r', label='d2ydx2')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('RBFs derivatives')
ax.legend()
fig.tight_layout()
Name = 'RBF_Derivatives.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

plt.show()

#%% RBF Solution procedure
Phi, Phi_x, Phi_xx = PHI_C4_X(x, x_b, c_r=0.1)
# Create the matrix L and the block A
L = Phi_xx + 2*Phi_x + Phi
A = 2*L.T@L
# Create the matrix Gamma_T
x_gamma = np.array([0, 1])
Phi_Gamma, Phi_x_Gamma, Phi_xx_Gamma = PHI_C4_X(x_gamma, x_b, c_r=0.1)
# Assembly the global system
Rows_H_1 = np.hstack((A, Phi_Gamma.T))
Rows_H_2 = np.hstack((Phi_Gamma, np.zeros((2, 2))))
A_star = np.vstack((Rows_H_1, Rows_H_2))
b_star = np.hstack((np.zeros(n_b), np.array([1, 3]))).T  
# Solve the system approximately using the Pseudo Inverse
x_sol = np.linalg.pinv(A_star, rcond=1e-15).dot(b_star)

# Get the norm of b for relative error:
r_error = np.linalg.norm(A_star.dot(x_sol)-b_star)
print(r_error)

# Get the weights
w = x_sol[:n_b]
# Assembly the solution
y_Sol_RBF = Phi.dot(w)
# Prepare the Analytic solution for comparison
y_Sol_Analyt = np.exp(-x) + (3*np.e-1) * x * np.exp(-x)

# This creates the figure
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.plot(x, y_Sol_RBF, 'b', label='RBF Solution')
ax.plot(x, y_Sol_Analyt, 'r', label='Analytic Solution')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Analytic vs RBF Solution')
fig.tight_layout()
Name = 'RBF_Approximation.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

plt.show()
