# -*- coding: utf-8 -*-
"""
Created on Jan 23 15:03:41 2024

@author: mendez
"""


import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("default", category=DeprecationWarning)
from tqdm import tqdm    
import shutil

import os
from matplotlib import cm

import imageio


#%% Preamble: customization of matplotlib
# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Rosenbrock function definition
def rosenbrock(X):  # The Function to Minimize  
    """Function to be minimized"""  # new        
    x = X[0]; y = X[1]
    C= 100*(y-x**2)**2+(1-x)**2.   
    return C 

#%% Function implementing the Simulated Annealing

def Simulated_A(func,X0,X_Bounds,n_iter,R,t_I,a=1):
    """
    Perform Simulated Annealing to minimize a given function.
    Parameters:
    func (callable): The objective function to minimize.
    X0 (array-like): Initial guess for the variables.
    X_Bounds (list of tuples): Bounds for each variable as (min, max).
    n_iter (int): Number of iterations to perform.
    R (float): Maximum change to apply to each variable in a step.
    t_I (float): Initial temperature.
    a (float, optional): Exponent in the cooling schedule. Default is 1.

    Returns:
    tuple: Arrays of the variable history, cost history, and the best solution found.
    """
    n_dim = len(X_Bounds)  # Number of dimensions
    X_next = np.zeros(n_dim)  # Next vector
    X_best = X0.copy()  # Best vector
    X_k = np.zeros((n_dim, n_iter + 1))  # Track positions
    Cost_k = np.zeros(n_iter + 1)  # Track costs

    # Initial evaluation
    X_k[:, 0] = X_best
    Cost_k[0] = func(X_best)
    F_best = Cost_k[0]

    print(f"Initial best solution: X = {X_best}, Cost = {F_best}")

    # Iterative search
    for k in range(1, n_iter + 1):
        # Generate a new candidate solution within bounds and exploration radius
        for j in range(n_dim):
            X_next[j] = np.random.uniform(
                max(X_k[j, k - 1] - R, X_Bounds[j][0]),
                min(X_k[j, k - 1] + R, X_Bounds[j][1]))
        # Evaluate the function at the new point
        F_next = func(X_next)
        # Update the best solution if the new one is better
        if F_next < F_best:
            X_best = X_next.copy()
            F_best = F_next
            print(f"Iteration {k}: Updated best solution to X = {X_best}, Cost = {F_best}")
        # Annealing process
        T = t_I / (k ** a + 1)
        diff = (F_next - F_best) if F_best != 0 else 0
        metropolis = np.exp(-diff / T) if diff > 0 else 1
        rand_change = np.random.rand() < metropolis
        # Accept new solution if better or based on Metropolis criterion
        if F_next < Cost_k[k - 1] or rand_change:
            X_k[:, k] = X_next
            Cost_k[k] = F_next
            if rand_change:
             print(f"Iteration {k}: Accepted random change, T={T:.2f}")
        else:
            X_k[:, k] = X_k[:, k - 1]
            Cost_k[k] = Cost_k[k - 1]

    return X_k, Cost_k, X_best

#%% Make animation of the RS History
def Anim(X_k,Cost_k,func,x_1m,x_1M,x_2m,x_2M,n_p,Name_Video):
   # This function makes an animation of the single search 
   # The X_k is the story of the optimizer.
   # func is the function optimized
   # x_1m,x_1M,x_2m,x_2M,n_p is the same input for the plot func
   # Name_Video is the name of the gif that will be exported
    
    # Temporary folder for storing images
    temp_folder = "Temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # Create a grid for contour plotting
    x = np.linspace(x_1m, x_1M, n_p)
    y = np.linspace(x_2m, x_2M, n_p)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    # Evaluate the objective function on the grid
    for i in range(n_p):
        for j in range(n_p):
            Z[i, j] = func(np.array([X[i, j], Y[i, j]]))

    # Get the minimum location in the grid
    min_loc = np.unravel_index(np.argmin(Z), Z.shape)

    # Generate images for each iteration
    n_iter = X_k.shape[1]
    for k in tqdm(range(n_iter),desc='Exporting Images'):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        # Plot the contour and optimization path
        ax1.contourf(X, Y, Z, cmap=cm.coolwarm, alpha=0.8)
        ax1.plot(X[min_loc], Y[min_loc], 'wo', markersize=5, label="Min")
        ax1.plot(X_k[0, :k + 1], X_k[1, :k + 1], 'ko-', label="Path")
        ax1.set_xlim([x_1m, x_1M])
        ax1.set_ylim([x_2m, x_2M])
        ax1.set_title(f"Iteration {k}")
        ax1.legend()

        # Plot the cost function history
        ax2.plot(range(k + 1), Cost_k[:k + 1], 'ro-')
        ax2.set_title("Cost History")
        ax2.set_xlabel("Iteration")
        ax2.set_ylabel("Cost")

        # Save the figure
        file_name = os.path.join(temp_folder, f"Step{k}.png")
        plt.savefig(file_name)
        plt.close(fig)

    # Create GIF
    images = []
    for k in range(n_iter):
        file_name = os.path.join(temp_folder, f"Step{k}.png")
        images.append(imageio.imread(file_name))

    # Save the GIF
    try:
        imageio.mimsave(Name_Video, images, duration=0.1)
        print(f"GIF successfully saved as {Name_Video}")
    except Exception as e:
        print(f"Error saving GIF: {e}")

    # Cleanup temporary folder
    if os.path.exists(temp_folder):
        shutil.rmtree(temp_folder)

#%% Compact Code
X_Bounds=[(-2,2),(-0.5,3)]  # Boundaries for x1 and x2
n_iter=200 # Number of iterations
X0=np.array([-2,0]) # Initial condition
t_I=rosenbrock(X0) # Initial temperature
R=0.2 # Width of the exploration
X_k,Cost_k,X_best =Simulated_A(rosenbrock,X0,X_Bounds,n_iter,R,t_I,1)

x_1m,x_1M,x_2m,x_2M,n_p=-2,2,-0.5,3,200
Anim(X_k,Cost_k,rosenbrock,x_1m,x_1M,x_2m,x_2M,n_p,'rosen_SA.gif')




