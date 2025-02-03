# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 10:35:23 2025

@author: mendez
"""

# This exercise is adapted from
# https://scikit-optimize.github.io/stable/auto_examples/bayesian-optimization.html

import numpy as np
np.random.seed(237)
import matplotlib.pyplot as plt
# This is a function to plot gaussian processes
from skopt.plots import plot_gaussian_process
from sklearn.metrics.pairwise import rbf_kernel
from scipy.linalg import cholesky, cho_solve
# Function for normalized gaussian
from scipy.stats import norm

import os



#%% Cost function definition
noise_level = 0.1
def f(x, noise_level=noise_level):
    return np.sin(5 * x[0]) * (1 - np.tanh(x[0] ** 2))\
           + np.random.randn() * noise_level

#%% Some home made functions 

# gp_fit: Fit the Gaussian Process
def gp_fit(x_s, y_s, l_c=0.3, sigma_y=0.1):
    # Compute the kernel matrix (training points only)
    K_ss = rbf_kernel(x_s, x_s, gamma=0.5 / l_c**2)
    # Cholesky decomposition
    L = cholesky(K_ss + sigma_y**2 * np.eye(len(x_s)), lower=True)
    # Solve for alpha
    alpha_v = cho_solve((L, True), y_s)
    return alpha_v, L

# gp_predict: Predictive mean and variance
def gp_predict(x, x_s, alpha_v, L, l_c=0.3, sigma_y=0.1):
    # Compute kernel matrices
    K = rbf_kernel(x, x, gamma=0.5 / l_c**2)
    K_s = rbf_kernel(x, x_s, gamma=0.5 / l_c**2)
    # Solve for K_alpha
    K_alpha = cho_solve((L, True), K_s.T)
    # Predictive mean and covariance
    mu_y = K_s @ alpha_v
    Sigma_yy = K - K_s @ K_alpha
    return mu_y, np.diag(Sigma_yy)


# Expected Improvement (EI) acquisition function
def expected_improvement(x, x_s, y_s, alpha_v, L, y_best, l_c=0.3, sigma_y=0.1, xi=0.01):
    # Reshape x for compatibility
    x = x.reshape(-1, 1)
    # Predictive mean and variance
    mu, sigma = gp_predict(x, x_s, alpha_v, L, l_c, sigma_y)
    # Normalize improvement
    with np.errstate(divide='warn'):
        Z = (y_best - mu - xi) / sigma
        ei = (y_best - mu - xi) * norm.cdf(Z) + sigma * norm.pdf(Z)
        ei[sigma == 0.0] = 0.0  # Handle cases where sigma is zero
    return ei


# Plot GP and acquisition function
def plot_gp_acquisition(x_grid,mu_y,uncertainty,EI,x_s,y_s,NAME_FIG):

    plt.figure(figsize=(14, 6))    
    # Left Plot: GP Regression
    plt.subplot(1, 2, 1)
    plt.fill_between(x_grid.ravel(),mu_y - uncertainty,  mu_y + uncertainty, 
        color="lightblue", alpha=0.5, label="95\% CI" )
    plt.plot(x_grid.ravel(), mu_y, "b-", label="GP Predictive Mean")
    plt.scatter(x_s, y_s, c="r", label="Observations")
    plt.axhline(np.min(y_s), color="g", linestyle="--", label="$y_{best}$")
    fx = [f(x_i, noise_level=0.0) for x_i in x_grid]
    plt.plot(x_grid.ravel(), fx, "r--", label="True (unknown)")
    plt.title("Gaussian Process Regression")
    plt.legend()

    # Right Plot: Acquisition Function
    plt.subplot(1, 2, 2)
    plt.plot(x_grid.ravel(), EI, "purple", label="EI Acquisition Function")
    plt.axvline(x_grid[np.argmax(EI)], color="red", linestyle="--", label="Proposed Point")
    plt.title("Acquisition Function")
    plt.legend()
    plt.show()
    plt.savefig(NAME_FIG, dpi=300, bbox_inches="tight")

#%% Run the BO and save the intermediate results
n_iter=10 # max number of iterations
n_ini=5 # initial number for the model definition
bounds=np.array([-2.0,2.0]) # bounds for the optimization
# Parameters of the GPR (very hard to guess!)
L_C=0.25
sigma_y=0.1
XI=0.01 # hyperparameter for the exploration

# define random initial samples:
x_init = np.random.uniform(bounds[0], bounds[1], size=(n_ini, 1))
y_init = np.array([f(x) for x in x_init])

# Initialize the set of points for the plotting
x_s,y_s= x_init,y_init
x_grid = np.linspace(bounds[0], bounds[1], 500).reshape(-1, 1)

# Prepare the convergence history
CONV=np.zeros(n_iter)


NAME_FOLDER='BO_RESULTS'
if not os.path.exists(NAME_FOLDER):
    os.makedirs(NAME_FOLDER)
    print(f"Folder '{NAME_FOLDER}' created.")

for i in range(n_iter): # BO loop
    # identify your current best solution
    y_best=np.min(y_s)
    CONV[i]=y_best
    # train the model on the current data available
    alpha_v, L = gp_fit(x_s, y_s, l_c=L_C, sigma_y=sigma_y) 
    # Show the current GPR on the grid (for plotting purposes only)
    mu_y, Sigma_y=gp_predict(x_grid,x_s, alpha_v, L, l_c=L_C, sigma_y=sigma_y )
    uncertainty = 1.96 * np.sqrt(Sigma_y)
    # Compute the EI on the current grid
    EI=expected_improvement(x_grid, x_s, y_s, alpha_v, L, y_best, l_c=L_C, sigma_y=sigma_y, xi=XI)
    # Plot and export the results
    NAME_FIG=f"{NAME_FOLDER}{os.sep}Iteration_{i}.png"
    plot_gp_acquisition(x_grid,mu_y,uncertainty,EI,x_s,y_s,NAME_FIG)    
    # Identify the next best location (sticking to the grid!)
    x_next=x_grid[np.argmax(EI)] # this step should be replaced by an optimization of EI      
    # Evaluate the cost function in the new proposal
    y_next = f(x_next)
    # Update the available dataset
    x_s = np.vstack((x_s, x_next))
    y_s = np.append(y_s, y_next)
    # Plot current GP and acquisition        
    print(f"Iteration {i + 1}: x_next = {x_next}, y_next = {y_next}")


fig, ax = plt.subplots(figsize=(5, 3)) 
plt.plot(CONV,'ko:')
ax.set_xlabel('Iterations',fontsize=16)
ax.set_ylabel('F',fontsize=16)  
plt.savefig('BO_Convergence.png',dpi=200,bbox_inches="tight")
plt.show()

plt.close('all')




