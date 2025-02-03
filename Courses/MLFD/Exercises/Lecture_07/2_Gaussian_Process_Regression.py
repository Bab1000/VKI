# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 19:45:45 2022

@author: mendez
"""


import numpy as np
import matplotlib.pyplot as plt

# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

#%% Define the problem to fit
x1 = np.linspace(0, 4.3, 200, endpoint=True)
x2 = np.linspace(4.8, 10, 200, endpoint=True)
x=np.concatenate((x1,x2))
# Create the deterministic part
y_clean= 3*x+(x/100)**3+4*np.sin(3/2*np.pi*x)
# Add (a seeded) stochastic part
np.random.seed(0)
y=y_clean+1*np.random.randn(len(x))
# Introduce some outliers in x=2 and x=8
G1=10*np.exp(-(x-2)**2/0.005)*np.random.randn(len(x))
G2=15*np.exp(-(x-8)**2/0.005)*np.random.randn(len(x))
y_final=y+G1+G2


#%% Normalize as usual

from sklearn.preprocessing import MinMaxScaler
# Scale the x 
scaler_X = MinMaxScaler(); 
scaler_X.fit_transform(x.reshape(-1,1))
x_prime=scaler_X.transform(x.reshape(-1,1)) # Scale

# Scale also the y 
scaler_Y = MinMaxScaler(); 
scaler_Y.fit_transform(y_final.reshape(-1,1))
y_prime=scaler_Y.transform(y_final.reshape(-1,1)) # Scale



#%% Use scikit learn's kernel functions
from sklearn.metrics.pairwise import rbf_kernel
def GP_regression(x,y,x_s,alpha=0.1,l_c=0.3):
    
    #1. Compute all portions of the covariance matrix
    K=rbf_kernel(x,x,gamma=0.5 / l_c**2)
    K_s=rbf_kernel(x,x_s,gamma=0.5 / l_c**2)
    K_ss=rbf_kernel(x_s,x_s,gamma=0.5 / l_c**2)
    
    #2. Prepare the required inversion
    K_inv=np.linalg.inv(K+alpha*np.eye(len(x)))
    
    #3. Compute mu and Sigma
    mu_y=K_s.T.dot(K_inv).dot(y)
    Sigma_y=K_ss-K_s.T@K_inv@K_s
    return mu_y, Sigma_y


#%% Plot GProcess with uncertainties
x_s=x_prime; y_s=y_prime
x_p=np.linspace(0,1,500).reshape(-1,1)
mu_y, Sigma_y=GP_regression(x_s,y_s,x_p,alpha=0.1,l_c=0.05)
uncertainty = 1.96 * np.sqrt(np.diag(Sigma_y))

fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_prime,y_prime,c='white',
            marker='o',edgecolor='black',
            s=10,label='Data')
plt.plot(x_p,mu_y,'r--')
plt.fill_between(x_p.ravel(), mu_y.ravel() + uncertainty.ravel(), \
                 mu_y.ravel() - uncertainty.ravel(), alpha=0.5)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)  
plt.savefig('GP_Regression.png',dpi=200)


#%% Implementation via scikit learn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
# Define the Process
rbf = RBF(length_scale=0.05)
gpr = GaussianProcessRegressor(kernel=rbf, alpha=0.01)
# Train the GP Regression
gpr.fit(x_s, y_s)
# Compute posterior mean and covariance; then uncertainties
mu_y, cov_s = gpr.predict(x_p, return_cov=True)
uncertainty = 1.96 * np.sqrt(np.diag(cov_s))

# Plot results as before
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_prime,y_prime,c='white',
            marker='o',edgecolor='black',
            s=10,label='Data')
plt.plot(x_p,mu_y,'r--')
plt.fill_between(x_p.ravel(), mu_y.ravel() + uncertainty.ravel(), \
                 mu_y.ravel() - uncertainty.ravel(), alpha=0.5)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)  
plt.savefig('GP_Regression_scikitlearn.png',dpi=200)












