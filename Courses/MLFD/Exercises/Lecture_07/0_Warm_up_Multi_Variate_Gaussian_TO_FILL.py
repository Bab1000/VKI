# -*- coding: utf-8 -*-
"""
Created on Sat May 25 19:35:28 2024

@author: mendez
"""

# We test the effect of marginalization and conditioning (slicing a Gaussian)


import numpy as np
import matplotlib.pyplot as plt # for plotting
from matplotlib import cm # colormaps for the plot

# Show the result of the fit 
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)


#%% Define mean and covariance
mu=np.array([2,1]) 
Sigma=np.array([[1,0.7],
               [0.7,2]])

# We re-use the same plotting function as in Exercise 4 of the pre-course

lam=np.max(np.linalg.eig(Sigma)[0])
n_g=100; # number of points in the grid

x_v=np.linspace(mu[0]-2*lam,mu[0]+2*lam,n_g)
y_v=np.linspace(mu[1]-2*lam,mu[1]+2*lam,n_g)
    
Xg1,Xg2=np.meshgrid(x_v,y_v)


def plot_2D_Gaussian(Xg1,Xg2,mu,Sigma):
  """
   Xg1,Xg2 define the grid over which the plot will be made.
   mu and Sigma are the mean and covariance matrix of the distribution
  """
  # note that the Gaussian is the exponential of a quadratic.
  # First we compute the shifted grid and Sigma's inverse:
  X_m_g1=Xg1-mu[0]; X_m_g2=Xg2-mu[1]    
  Sigma_inv=np.linalg.inv(Sigma)
  Det_sigma=np.linalg.det(Sigma)  
  Scaling=1/np.pi*(1/(np.sqrt(2)*Det_sigma))
  
  # then we loop over the grid (not efficient, but ok!)
  N=np.zeros_like(Xg1)
  for i in range(np.shape(Xg1)[0]):
   for j in range(np.shape(Xg1)[1]):
    x_m=np.array([X_m_g1[i,j],X_m_g2[i,j]])   
    N[i,j]=Scaling*np.exp(-x_m.T.dot(Sigma_inv).dot(x_m))
  
  plt.figure()
  plt.contourf(Xg1,Xg2,N,cmap=cm.plasma)



# Plot it ! 
plt.figure(figsize=(4, 4))
plot_2D_Gaussian(Xg1,Xg2,mu,Sigma)  
  

#%% Step 1: plot the marginal (fill me!)
mu_1=FILL ME
sigma_1=FILL ME

# Take a range of 3 sigmas
x = np.linspace(mu_1 - 3*sigma_1, mu_1 + 3*sigma_1, 100)

# This is the resulting Gaussian
y = (1 / (sigma_1 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu_1) / sigma_1)**2)
    
# Creating the plot
plt.figure(figsize=(4, 3))
plt.plot(x, y)
plt.title('$p(x_1)$')
plt.xlabel('$x_1$')
plt.ylabel('Marginal Probability Density')
plt.tight_layout()
Name_FIG='Distribution_1'
plt.savefig(Name_FIG)
plt.show()


#%% Step 2: plot the marginal (fill me!)
mu_1=FILL ME
sigma_1= FILL ME

# Take a range of 3 sigmas
x = np.linspace(mu_1 - 3*sigma_1, mu_1 + 3*sigma_1, 100)

# This is the resulting Gaussian
y = (1 / (sigma_1 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu_1) / sigma_1)**2)
    
# Creating the plot
plt.figure(figsize=(4, 3))
plt.plot(x, y)
plt.title('$p(x_1)$')
plt.xlabel('$x_1$')
plt.ylabel('Marginal Probability Density')
plt.tight_layout()
Name_FIG='Distribution_2'
plt.savefig(Name_FIG)
plt.show()
















