# -*- coding: utf-8 -*-
"""
Forward Propagation of a Linearized Pendulum with Uncertainty
Created on Thu Jan 16 13:48:28 2025
Author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.linalg import expm  # For computing matrix exponentials
from scipy.integrate import odeint



# Customize plot settings for LaTeX and larger fonts
plt.rcParams.update({
    'text.usetex': True,
    'font.size': 16,
    'font.family': 'serif'
})

# =============================================================================
# Load and Preprocess the Measurement Data
# =============================================================================

# Measurement file number
MEAS = 1
file_name = f'data_pendulum/measurement_{MEAS}.dat'

# Load measurement data
data = pd.read_csv(file_name, delimiter='\t', skiprows=1)
t_s = data.values[:, 0]  # Time vector
s1_s = data.values[:, 1]  # Measured angle (theta)
s2_s = data.values[:, 2]  # Measured angular velocity (theta_dot)
s_s = np.array([s1_s, s2_s])  # Combined state vector

# Extract time step and number of points
n_t = len(t_s)
dt = t_s[2] - t_s[1]

# =============================================================================
# Define System Parameters and Matrices
# =============================================================================

# Linearized system parameters
omega_n = 5.73  # Natural frequency
mu = 0.18       # Damping coefficient

# Continuous-time system matrix (linearized pendulum dynamics)
A = np.array([
    [0, 1],
    [-omega_n**2, -mu]
])

# Compute discrete-time state transition matrix using exact discretization
A_d = expm(A * dt)  # Discrete-time state transition matrix
# test the Euler approximation also: (it does not work!)
# A_d=np.eye(2)+A*dt

#%% =============================================================================
# Initialize State and Covariance Matrices
# =============================================================================
# Initial state from the measurement data
x_0 = np.array([s1_s[0], s2_s[0]])  # Initial angle and angular velocity
# Initialize state and covariance matrices for propagation
s_p = np.zeros_like(s_s)  # Predicted state over time
Sigma_p = np.zeros((2, 2, n_t))  # Predicted covariance matrices over time
# Initial conditions
s_p[:, 0] = x_0
# Noise covariance matrix (process noise)
Q = np.array([[0.00002, 0],[0, 0.0002]])
Sigma_p[:, :, 0] = Q  # Initial uncertainty guess

# Measurement noise (try larger values as well)
R = np.array([[0.02, 0],[0, 0.02]])


# Store uncertainties for plotting (95% CI)
delta_s1 = np.zeros(n_t)  # Uncertainty in angle
delta_s2 = np.zeros(n_t)  # Uncertainty in angular velocity
delta_s1[0] = 1.96 * np.sqrt(Sigma_p[0, 0, 0])
delta_s2[0] = 1.96 * np.sqrt(Sigma_p[1, 1, 0])

#%% =============================================================================
# Forward Propagation with Uncertainty
# =============================================================================
for k in range(n_t - 1):
    # Prediction step
    s_pred = A_d @ s_p[:, k]  # Predict the next state
    # Predict the next covariancee
    Sigma_pred = A_d @ Sigma_p[:, :, k] @ A_d.T + Q  
    # Kalman gain
    K = Sigma_pred @ np.linalg.inv(Sigma_pred  + R)
    # Update step
    z_kp1 = s_s[:, k + 1]  # Measurement at step k+1
    s_p[:, k + 1] = s_pred + K @ (z_kp1 - s_pred)  # Update state estimate
    Sigma_p[:, :, k + 1] = (np.eye(2) - K) @ Sigma_pred  # Update covariance
    # Compute uncertainties (95% CI)
    delta_s1[k + 1] = 1.96 * np.sqrt(Sigma_p[0, 0, k + 1])
    delta_s2[k + 1] = 1.96 * np.sqrt(Sigma_p[1, 1, k + 1])
# =============================================================================
# Visualization of Results
# =============================================================================

# we compute also the true solution:
degs=np.array([+35,-28,+15,-45,36,-22,18,42])
s0_s=degs*np.pi/180
mu_s=np.array([0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18])

mu=mu_s[MEAS]      
s0=[s0_s[MEAS],0]
    
# We define the function of the system    
def f_ODEINT(s, t,mu=0.5,omega_n=5):
    # this is the function f in the slides
    s_1_dot = s[1]
    s_2_dot = -mu*s[1]-omega_n**2*np.sin(s[0])
    return [s_1_dot, s_2_dot]


Y = odeint(f_ODEINT, s0, t_s,args=(mu,omega_n))
s1_f=Y[:,0]; s2_f=Y[:,1] # Get the solution
    




# Plot angle (theta) with uncertainty
fig, ax = plt.subplots(figsize=(8, 5))
plt.plot(t_s, s_p[0, :], label="Predicted $s_1$ (Angle)")
plt.plot(t_s, s_s[0, :], 'r--', label="Measured $s_1$")
plt.plot(t_s,s1_f,'k',label="Truth")
plt.fill_between(t_s,
                 s_p[0, :] - delta_s1,
                 s_p[0, :] + delta_s1,
                 color='blue', alpha=0.2, label="95\% CI for $s_1$")
plt.ylim([1.5 * np.min(s_p[0, :]), 1.5 * np.max(s_p[0, :])])
plt.xlabel('$t$', fontsize=18)
plt.ylabel('$\\theta(t)$', fontsize=18)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
plt.savefig(f'Demo_Theta_Evolution_FULL_{MEAS}.png', dpi=300)

# Plot angular velocity (theta_dot) with uncertainty
fig, ax = plt.subplots(figsize=(8, 5))
plt.plot(t_s, s_p[1, :], label="Predicted $s_2$ (Angular Velocity)")
plt.plot(t_s, s_s[1, :], 'r--', label="Measured $s_2$")
plt.plot(t_s,s2_f,'k',label="Truth")
plt.fill_between(t_s,
                 s_p[1, :] - delta_s2,
                 s_p[1, :] + delta_s2,
                 color='blue', alpha=0.2, label="95\% CI for $s_2$")
plt.ylim([1.5 * np.min(s_p[1, :]), 1.5 * np.max(s_p[1, :])])
plt.xlabel('$t$', fontsize=18)
plt.ylabel('$\dot{\\theta}(t)$', fontsize=18)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
plt.savefig(f'Demo_Theta_Dot_Evolution_FULL_{MEAS}.png', dpi=300)
