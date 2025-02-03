# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:15:50 2025

@author: mendez, ratz
"""

# Utilities functions for setting up the data and plotting
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import qmc
from mpl_toolkits.axes_grid1 import make_axes_locatable

# This imports a reduced version of the SPICY class which is based on the
# original publication:
# Sperotto, Pieraccini, Mendez (2022) https://arxiv.org/abs/2112.12752
# An open source implementation can be found here:
# Sperotto, Ratz, Mendez (2024) https://joss.theoj.org/papers/10.21105/joss.05749#
from spicy_class import spicy

Fol_Out = '04_2D_Regression_Exercise'
if not os.path.exists(Fol_Out):
    os.makedirs(Fol_Out)

# Configuration for plots  
fontsize = 12
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['legend.fontsize'] = fontsize-3
plt.rcParams['font.size'] = fontsize

# Fix the random seed
np.random.seed(42)

# Set the domain boundaries
x_min, x_max = -1, 1
y_min, y_max = 0, 1

# Set the properties of the corner flow
A = 1
n = 4/3
beta = np.pi/n

# Get the number of points in the domain. The actual number of points is smaller
# since we are looking at a corner flow. For n = 4/3, the actual number of points
# is approximately n_p_full*0.75
n_p_full = 400
# We use Halton points since they are semirandomly spread in space
sampler = qmc.Halton(d=2, scramble=False)
sample = sampler.random(n=n_p_full).T
# We rescale them to the domain boundaries
X_train = sample[0] * (x_max-x_min) + x_min
Y_train = sample[1] * (y_max-y_min) + y_min

# Here, we filter the points which are not inside the domain of the corner flow
valid_lower = Y_train >= 0
valid_incli = Y_train >= X_train*np.tan(beta)*0.999 # needed because of numerical precision
valid = valid_lower * valid_incli
X_train = X_train[valid]
Y_train = Y_train[valid]
n_p = X_train.shape[0]

# This computes the analytical solution based on potential flow
r = np.sqrt(X_train**2 + Y_train**2)
theta = np.arctan2(Y_train, X_train)
# Get the radial and tangential velocity
U_r = A * n * r**(n-1) * np.cos(n * theta)
U_t = -A * n * r**(n-1) * np.sin(n * theta)
# Get the x- and y-velocity
U_train = U_r*np.cos(theta) - U_t*np.sin(theta)
V_train = U_r*np.sin(theta) + U_t*np.cos(theta)

# We add uniform noise to the data to simulate measurement noise. This value
# is very extreme, even for measurements
q = 0.3
U_noise = U_train * (1 + q*np.random.uniform(-1, 1, n_p))
V_noise = V_train * (1 + q*np.random.uniform(-1, 1, n_p))

# Plot the noisy, raw data as a quiver plot
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
ax.quiver(X_train, Y_train, U_noise, V_noise, pivot='middle', label=r'Training data $\bm{X}$')
ax.plot([0, 1], [0, 0], c='r', lw=3, label='Walls')
ax.plot([0, -1], [0, 1], c='r', lw=3)
ax.set_aspect(1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Training_data.png')
plt.show()


#%% Here, we start to set up the regression with the spicy class

# These are Dirichlet boundary conditions on the horizontal wall
angle_hor = 0
n_Dir_hor = 50
X_Dir_hor = np.linspace(0, x_max, n_Dir_hor)
Y_Dir_hor = np.zeros(n_Dir_hor)
V_Dir_hor = np.zeros(n_Dir_hor)
n_x_hor = np.ones(n_Dir_hor) * np.sin(angle_hor)
n_y_hor = np.ones(n_Dir_hor) * np.cos(angle_hor)

# These are Dirichlet boundary conditions on the inclined wall
angle_incl = beta
n_Dir_incl = 50
Y_Dir_incl = np.linspace(y_min, y_max, n_Dir_incl)
X_Dir_incl = Y_Dir_incl / np.tan(beta)
U_Dir_incl = np.zeros(n_Dir_incl)
n_x_incl = np.ones(n_Dir_incl) * np.sin(angle_incl)
n_y_incl = -np.ones(n_Dir_incl) * np.cos(angle_incl)


# Here, we set up the curl-free constraints (this is a regular grid that we 
# use for the out of sample results)
n_pred = 31
x_pred = np.linspace(x_min, x_max, n_pred*2-1)
y_pred = np.linspace(y_min, y_max, n_pred)
X_pred, Y_pred = np.meshgrid(x_pred, y_pred)
X_pred = X_pred.ravel()
Y_pred = Y_pred.ravel()
# Filter them to be inside the domain
valid_lower = Y_pred >= 0
valid_incli = Y_pred >= X_pred*np.tan(beta)*0.999
valid = valid_lower * valid_incli
X_pred = X_pred[valid]
Y_pred = Y_pred[valid]

# Set up the lists for the constraint coordinates. Curl and Divergence are 
# constrained in the training and prediction data. The Walls are constrained
# with Dirichlet constraints. From the horizontal wall, we remove the first point
# because it is already included in the inclined wall. Repeating the constraint
# leads to horrible conditioning of the basis!
CURL = [np.concatenate((X_train, X_pred)), np.concatenate((Y_train, Y_pred))]
DIV = [np.concatenate((X_train, X_pred)), np.concatenate((Y_train, Y_pred))]
DIR_hor = [X_Dir_hor[1:], Y_Dir_hor[1:], n_x_hor[1:], n_y_hor[1:], V_Dir_hor[1:]]
DIR_incl = [X_Dir_incl, Y_Dir_incl, n_x_incl, n_y_incl, U_Dir_incl]


# Initialize the spicy object with the data
SP = spicy(data=[U_noise, V_noise], grid_point=[X_train, Y_train])
# Perform the clustering to set the collocation points
SP.clustering(n_K=[1.1], Areas=[[]], r_mM=[0.03, 0.3], eps_l=0.88)
# Set the boundary conditions
SP.vector_constraints(DIR=[DIR_hor, DIR_incl], CURL=CURL, DIV=DIV, extra_RBF=True)

#%% This code block visualizes the clustering result as well as the constraints

fig, ax = plt.subplots(figsize = (5, 3), dpi=200)
# Get the collocation points from the data
X_Plot = SP.X_C[np.argwhere(SP.Clust_list==0)]
Y_Plot = SP.Y_C[np.argwhere(SP.Clust_list==0)]
d_K_Plot = SP.d_k[np.argwhere(SP.Clust_list==0)]
# Plot collocation circles from the data
for i in range(0,len(X_Plot),1):
    circle1 = plt.Circle((X_Plot[i], Y_Plot[i]), d_K_Plot[i]/2, lw=2,
                          fill=True, facecolor='g', edgecolor='black', alpha=0.05)
    ax.add_artist(circle1)
# Plot the data
ax.scatter(SP.X_G, SP.Y_G, c=np.sqrt(SP.u**2 + SP.v**2), edgecolor='black')
ax.set_aspect(1)
ax.set_title("RBF Collocation for data points")
ax.set_xlabel('x')
ax.set_ylabel('y')
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Collocation_points_data.png')


fig, ax = plt.subplots(figsize = (5, 3), dpi=200)
# Get the collocation points from the constraints
X_Plot = SP.X_C[np.argwhere(SP.Clust_list==1)]
Y_Plot = SP.Y_C[np.argwhere(SP.Clust_list==1)]
d_K_Plot = SP.d_k[np.argwhere(SP.Clust_list==1)]
# Plot collocation circles from the constraints
for i in range(0,len(X_Plot),1):
    circle1 = plt.Circle((X_Plot[i], Y_Plot[i]), d_K_Plot[i]/2, lw=0.5,
                          fill=True, facecolor='g', edgecolor='k', alpha=0.1)
    ax.add_artist(circle1)
# Plot the constraints
ax.plot(SP.X_D[0], SP.Y_D[0], 'ro', label='Dirichlet, bottom')
ax.plot(SP.X_D[1], SP.Y_D[1], 'r*', label='Dirichlet, inclined')
ax.plot(SP.X_Div, SP.Y_Div, 'bd', label='Divergence')
ax.plot(SP.X_Curl, SP.Y_Curl, 'bv', label='Curl')
ax.set_aspect(1)
ax.set_title("RBF Collocation for constraints")
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()
# Save the result
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Collocation_points_constraints.png')

plt.show()

#%% Now, we can assemble the linear system using a penalty on the divergence
# as additional regularization
SP.Assembly_Regression(alpha_div=1)
# Solve the system with a fixed conditioning number
SP.Solve(K_cond=1e6)
# Get the solution and derivatives in the training points and print their errors
U_train_RBF, V_train_RBF = SP.Get_Sol([X_train, Y_train])


dUdX_train_RBF, dUdY_train_RBF, dVdX_train_RBF, dVdY_train_RBF = SP.Get_first_Derivatives([X_train, Y_train])

error_u = np.linalg.norm(U_train_RBF - U_train) / np.linalg.norm(U_train)
error_v = np.linalg.norm(V_train_RBF - V_train) / np.linalg.norm(V_train)

print('In-sample velocity error in u: {0:.3f}%'.format(error_u*100))
print('In-sample velocity error in v: {0:.3f}%'.format(error_v*100))
print('In-sample divergence error: {0:.3f}'.format(np.sum(np.abs(dUdX_train_RBF + dVdY_train_RBF))))
print('In-sample curl error: {0:.3f}'.format(np.sum(np.abs(dUdY_train_RBF - dVdX_train_RBF))))

# Get the solution and derivatives in the prediction points and print their errors
# First, generate the ground truth
r_pred = np.sqrt(X_pred**2 + Y_pred**2)
theta_pred = np.arctan2(Y_pred, X_pred)
U_r_pred = A * n * r_pred**(n-1) * np.cos(n * theta_pred)
U_t_pred = -A * n * r_pred**(n-1) * np.sin(n * theta_pred)
U_pred = U_r_pred*np.cos(theta_pred) - U_t_pred*np.sin(theta_pred)
V_pred = U_r_pred*np.sin(theta_pred) + U_t_pred*np.cos(theta_pred)

U_pred_RBF, V_pred_RBF = SP.Get_Sol([X_pred, Y_pred])
dUdX_pred_RBF, dUdY_pred_RBF, dVdX_pred_RBF, dVdY_pred_RBF = SP.Get_first_Derivatives([X_pred, Y_pred])

error_u_pred = np.linalg.norm(U_pred - U_pred_RBF) / np.linalg.norm(U_pred)
error_v_pred = np.linalg.norm(V_pred - V_pred_RBF) / np.linalg.norm(V_pred)

print('Out-of-sample velocity error in u: {0:.3f}%'.format(error_u_pred*100))
print('Out-of-sample velocity error in v: {0:.3f}%'.format(error_v_pred*100))
print('Out-of-sample divergence error: {0:.3f}'.format(np.sum(np.abs(dUdX_pred_RBF+dVdY_pred_RBF))))
print('Out-of-sample curl error: {0:.3f}'.format(np.sum(np.abs(dUdY_pred_RBF-dVdX_pred_RBF))))

#%% Now, we can look at the detailed results in plots

# Plot the quiver plot of the training and prediction data

fig, ax = plt.subplots(figsize = (5, 3), dpi=200)
ax.set_aspect(1)
ax.quiver(X_train, Y_train, U_train, V_train, pivot='middle')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title(r'Training data')
ax.set_aspect(1)
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_quiver_plot_training.png')


fig, ax = plt.subplots(figsize = (5, 3), dpi=200)
ax.set_aspect(1)
ax.quiver(X_pred, Y_pred, U_pred, V_pred, pivot='middle')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title(r'Prediction data')
ax.set_aspect(1)
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_quiver_plot_prediction.png')
plt.show()

#%%


# Plot a contour of the difference in the training data for the velocities
fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_train, Y_train, c=np.abs(U_train - U_train_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.3)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Velocity error $u$: Training data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_velocity_error_u_x_training.png')
plt.show()



fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_train, Y_train, c=np.abs(V_train - V_train_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Velocity error $v$: Training data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_velocity_error_v_x_training.png')
plt.show()


# Plot a contour of the difference in the training data for curl and divergence
fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_train, Y_train, c=np.abs(dUdX_train_RBF + dVdY_train_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Divergence error: Training data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_div_x_training.png')
plt.show()

fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_train, Y_train, c=np.abs(dUdY_train_RBF - dVdX_train_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Curl error: Training data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_curl_x_training.png')
plt.show()


# Plot a contour of the difference in the prediction data for the velocities
fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_pred, Y_pred, c=np.abs(U_pred_RBF - U_pred))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Velocity error $u$: Prediction data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_velocity_error_u_x_prediction.png')
plt.show()

fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_pred, Y_pred, c=np.abs(V_pred_RBF - V_pred))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Velocity error $v$: Prediction data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_velocity_error_v_x_prediction.png')
plt.show()


# Plot a contour of the difference in the prediction data for curl and divergence
fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_pred, Y_pred, c=np.abs(dUdX_pred_RBF+dVdY_pred_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Divergence error: Prediction data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_div_x_prediction.png')
plt.show()

fig, ax = plt.subplots(figsize = (5, 3.5), dpi=200)
ax.set_aspect(1)
cont = ax.scatter(X_pred, Y_pred, c=np.abs(dUdY_pred_RBF-dVdX_pred_RBF))
cax = make_axes_locatable(ax).append_axes('top', size='5%', pad=0.4)
cbar = fig.colorbar(cont, cax=cax, orientation='horizontal')
cbar.set_label(r'Curl error: Prediction data', labelpad=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.tight_layout()
fig.savefig(Fol_Out + os.sep + 'Result_curl_x_prediction.png')
plt.show()